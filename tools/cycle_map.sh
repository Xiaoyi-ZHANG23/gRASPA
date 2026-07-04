#!/usr/bin/env bash
# cycle_map.sh — map field convergence (Tanimoto) vs production cycles for the
# unconverged pilot MOFs, to build a cost-per-converged-field model vs loading.
# Runs each MOF at several PROD cycle counts x a seed pair, across the 8 H100s.
set -uo pipefail
NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/24.5
export LD_LIBRARY_PATH="$NVHPC/compilers/lib:$NVHPC/cuda/lib64:$NVHPC/math_libs/lib64:${LD_LIBRARY_PATH:-}"
BIN=/home/nvidia/gRASPA/src_clean/nvc_main.x
INPUT_DIR=/home/nvidia/Smoke_Test/INPUT
CIF_DIR=/home/nvidia/Smoke_Test/target
TOOLS=/home/nvidia/gRASPA/tools
OUT=/data/pilot/cyclemap
COSTCSV=/data/pilot/pilot_costs.csv
NGPU=8; INIT=50000; SEEDS="1 999"

# 5 unconverged MOFs (name + unit cells), loading noted for reference
MOFS=(
  "qmof-d475ee7 1 1 1"    # load 154
  "qmof-b334eed 3 3 3"    # load 192
  "qmof-b4945bb 2 2 2"    # load 245
  "qmof-b9c6168 3 3 2"    # load 265
  "qmof-7e3c643 2 2 2"    # load 327
)
PRODS="4000000 2000000 1000000"   # longest first for packing

# build job list (longest prod first)
JOBS=()
for prod in $PRODS; do for m in "${MOFS[@]}"; do for s in $SEEDS; do
  JOBS+=("$m|$prod|$s"); done; done; done

mkdir -p "$OUT"
declare -a SLOT_PID
run_one () { # $1 gpu ; $2 "name nx ny nz" ; $3 prod ; $4 seed
  local g="$1" mof="$2" prod="$3" seed="$4"
  local name nx ny nz; read -r name nx ny nz <<<"$mof"
  local R="$OUT/${name}/p${prod}_s${seed}"; rm -rf "$R"; mkdir -p "$R"
  cp "$INPUT_DIR"/simulation.input "$INPUT_DIR"/*.def "$R"/ 2>/dev/null
  cp "$CIF_DIR/$name.cif" "$R"/ || { echo "MISSING $name"; return 1; }
  sed -i "s/^FrameworkName .*/FrameworkName $name/;
          s/^UnitCells .*/UnitCells 0 $nx $ny $nz/;
          s/^RandomSeed .*/RandomSeed  $seed/;
          s/^NumberOfInitializationCycles .*/NumberOfInitializationCycles $INIT/;
          s/^NumberOfProductionCycles .*/NumberOfProductionCycles     $prod/" "$R"/simulation.input
  ( cd "$R" && CUDA_VISIBLE_DEVICES="$g" "$BIN" > output.txt 2> run_stderr.txt )
  # loading/block-averages default to stdout since upstream 62b25fc -> parse output.txt
  python3 "$TOOLS/field_cost.py" --run-dir "$R" --log output.txt \
          --replicate "$seed" --gpu-id "$g" --gpus 1 --csv "$COSTCSV" \
          --notes "cyclemap $name p$prod s$seed" >/dev/null 2>&1
  echo "  done: $name p$prod s$seed (gpu$g)"
}
dispatch () {
  while true; do
    for g in $(seq 0 $((NGPU-1))); do
      local pid=${SLOT_PID[$g]:-}
      if [[ -z "$pid" ]] || ! kill -0 "$pid" 2>/dev/null; then
        run_one "$g" "$@" & SLOT_PID[$g]=$!; return; fi
    done; sleep 3
  done
}
echo "=== cycle_map: ${#JOBS[@]} runs (5 MOFs x {1M,2M,4M} x seeds 1/999) ==="
for job in "${JOBS[@]}"; do IFS='|' read -r mof prod seed <<<"$job"; dispatch "$mof" "$prod" "$seed"; done
wait

echo "=== convergence vs cycles (continuous Tanimoto, .cube) ==="
CONV="$OUT/cyclemap_convergence.csv"
echo "framework,prod_cycles,tanimoto,converged" > "$CONV"
for m in "${MOFS[@]}"; do read -r name _ <<<"$m"
  for prod in 1000000 2000000 4000000; do
    line=$(python3 "$TOOLS/field_convergence.py" --run-dirs "$OUT/$name/p${prod}_s1" "$OUT/$name/p${prod}_s999" --threshold 0.75 2>/dev/null)
    tan=$(grep -oE 'Tanimoto=[0-9.]+' <<<"$line"|head -1|cut -d= -f2)
    ok=$(grep -qE 'CONVERGED \(PASS\)' <<<"$line" && echo yes || echo no)
    echo "$name,$prod,$tan,$ok" | tee -a "$CONV"
  done
done
echo "=== done -> $CONV ==="
