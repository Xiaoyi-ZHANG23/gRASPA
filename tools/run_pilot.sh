#!/usr/bin/env bash
# run_pilot.sh — single-component Xe density-field pilot across the local 8x H100 node.
#
# No SLURM here: this dispatches runs directly across GPUs 0..NGPU-1 with a slot
# semaphore (one sim per GPU at a time). For each MOF x seed it stages a run dir,
# substitutes FrameworkName / UnitCells / RandomSeed / cycle counts into a copy of
# the Xe simulation.input, runs gRASPA, logs per-field cost, and finally gates field
# convergence on each MOF's replicate pair (using the triclinic-correct .cube).
#
# USAGE:
#   tools/run_pilot.sh --dry-run          # show the plan + GPU-hour estimate, run nothing
#   tools/run_pilot.sh                     # run the pilot
#   PROD=500000 SEEDS="1 999" tools/run_pilot.sh --dry-run
#
# Everything below is overridable via environment variables.
set -uo pipefail

# ---- config ----------------------------------------------------------------
NVHPC="${NVHPC:-/opt/nvidia/hpc_sdk/Linux_x86_64/24.5}"
BIN="${BIN:-/home/nvidia/gRASPA/src_clean/nvc_main.x}"
INPUT_DIR="${INPUT_DIR:-/home/nvidia/Smoke_Test/INPUT}"          # Xe simulation.input + *.def
CIF_DIR="${CIF_DIR:-/home/nvidia/Smoke_Test/target}"            # qmof-*.cif live here
MOF_LIST="${MOF_LIST:-/home/nvidia/Smoke_Test/pilot_mofs.txt}"
OUTROOT="${OUTROOT:-/data/pilot/xe_pilot}"                      # durable storage
COSTCSV="${COSTCSV:-/data/pilot/pilot_costs.csv}"
CONVCSV="${CONVCSV:-$OUTROOT/pilot_convergence.csv}"
TOOLS="${TOOLS:-/home/nvidia/gRASPA/tools}"
NGPU="${NGPU:-8}"
SEEDS="${SEEDS:-1 999}"          # DISTINCT NONZERO seeds (glibc: 0==1, never pair them)
INIT="${INIT:-50000}"
PROD="${PROD:-500000}"          # production cycles (set from the convergence ladder)
USCYCLE="${USCYCLE:-233}"       # measured microseconds/cycle for the estimate (UseMaxStep yes)

export LD_LIBRARY_PATH="$NVHPC/compilers/lib:$NVHPC/cuda/lib64:$NVHPC/math_libs/lib64:${LD_LIBRARY_PATH:-}"

# ---- read MOF list ---------------------------------------------------------
mapfile -t ROWS < <(grep -vE '^\s*#|^\s*$' "$MOF_LIST")
JOBS=()
for row in "${ROWS[@]}"; do
  name=$(awk '{print $1}' <<<"$row"); nx=$(awk '{print $2}' <<<"$row")
  ny=$(awk '{print $3}' <<<"$row");  nz=$(awk '{print $4}' <<<"$row")
  for s in $SEEDS; do JOBS+=("$name|$nx $ny $nz|$s"); done
done

# ---- estimate --------------------------------------------------------------
NMOF=${#ROWS[@]}; NJOB=${#JOBS[@]}; TOTCYC=$((INIT+PROD))
sec_per_run=$(python3 -c "print($TOTCYC*$USCYCLE/1e6)")
gpuh=$(python3 -c "print(round($NJOB*$sec_per_run/3600,3))")
wall=$(python3 -c "import math;print(round(math.ceil($NJOB/$NGPU)*$sec_per_run/60,1))")
echo "=== Xe pilot plan ==="
echo "MOFs: $NMOF   seeds: [$SEEDS]   runs: $NJOB   init/prod: $INIT/$PROD (${TOTCYC} cyc)"
echo "est per-run: ${sec_per_run}s @ ${USCYCLE}us/cyc   total: ${gpuh} GPU-hours   ~${wall} min wall on ${NGPU} GPUs"
echo "outputs -> $OUTROOT (durable)"
if [[ "${1:-}" == "--dry-run" ]]; then
  printf '  %s\n' "${JOBS[@]/|/ }"; exit 0
fi

mkdir -p "$OUTROOT"
# ---- GPU-slot semaphore ----------------------------------------------------
declare -a SLOT_PID
run_one () {  # $1=gpu $2=name $3="nx ny nz" $4=seed
  local g="$1" name="$2" uc="$3" seed="$4"
  local R="$OUTROOT/$name/s$seed"; rm -rf "$R"; mkdir -p "$R"
  cp "$INPUT_DIR"/simulation.input "$INPUT_DIR"/*.def "$R"/ 2>/dev/null
  cp "$CIF_DIR/$name.cif" "$R"/ || { echo "MISSING CIF $name"; return 1; }
  sed -i "s/^FrameworkName .*/FrameworkName $name/;
          s/^UnitCells .*/UnitCells 0 $uc/;
          s/^RandomSeed .*/RandomSeed  $seed/;
          s/^NumberOfInitializationCycles .*/NumberOfInitializationCycles $INIT/;
          s/^NumberOfProductionCycles .*/NumberOfProductionCycles     $PROD/" "$R"/simulation.input
  ( cd "$R" && CUDA_VISIBLE_DEVICES="$g" "$BIN" > output.txt 2> run_stderr.txt )
  python3 "$TOOLS/field_cost.py" --run-dir "$R" --log run_stderr.txt \
          --replicate "$seed" --gpu-id "$g" --gpus 1 --csv "$COSTCSV" \
          --notes "xe_pilot $name s$seed" >/dev/null 2>&1
  echo "  done: $name s$seed (gpu$g)"
}
dispatch () {  # find a free GPU slot, launch, block only if all busy
  while true; do
    for g in $(seq 0 $((NGPU-1))); do
      local pid=${SLOT_PID[$g]:-}
      if [[ -z "$pid" ]] || ! kill -0 "$pid" 2>/dev/null; then
        run_one "$g" "$@" & SLOT_PID[$g]=$!; return
      fi
    done; sleep 2
  done
}
echo "=== launching $NJOB runs ==="
for job in "${JOBS[@]}"; do IFS='|' read -r name uc seed <<<"$job"; dispatch "$name" "$uc" "$seed"; done
wait
echo "=== all runs done; gating convergence (.cube, Tanimoto>=0.75) ==="
echo "framework,seed_a,seed_b,tanimoto,binary_tanimoto,js_bits,converged" > "$CONVCSV"
read -r -a SARR <<<"$SEEDS"
for name in $(awk '{print $1}' <<<"$(grep -vE '^\s*#|^\s*$' "$MOF_LIST")"); do
  a="$OUTROOT/$name/s${SARR[0]}"; b="$OUTROOT/$name/s${SARR[1]}"
  line=$(python3 "$TOOLS/field_convergence.py" --run-dirs "$a" "$b" --threshold 0.75 --binary-cutoff 0.05 2>/dev/null)
  tan=$(grep -oE 'Tanimoto=[0-9.]+' <<<"$line"|head -1|cut -d= -f2)
  bt=$(grep -oE 'binT=[0-9.]+' <<<"$line"|head -1|cut -d= -f2)
  js=$(grep -oE 'JS=[0-9.]+' <<<"$line"|head -1|cut -d= -f2)
  ok=$(grep -qE 'CONVERGED \(PASS\)' <<<"$line" && echo yes || echo no)
  echo "$name,${SARR[0]},${SARR[1]},$tan,$bt,$js,$ok" | tee -a "$CONVCSV"
done
echo "=== convergence -> $CONVCSV ; costs -> $COSTCSV ==="
