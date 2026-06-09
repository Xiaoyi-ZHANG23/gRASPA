#!/usr/bin/env python3
"""
score.py — Parse gRASPA output and compare against a reference baseline.

VENDORED, unmodified, from AutoJIT-gRASPA (https://github.com/Zhaoli2042/AutoJIT-gRASPA),
scripts/score.py. Bundled here so the gRASPA debugging workflow (debugging/REPRODUCE.md) is
reproducible from this repo alone, without a separate AutoJIT clone. See debugging/README.md for
the comparison-class vs absolute-class check distinction (it matters when reading the verdict).

Correctness checks (all must pass):
  1. Move statistics (exact match — same RNG = same moves)
  2. PseudoAtom counts (exact match — same RNG = same final state)
  3. Final energies (match within tolerance vs reference)
  4. Energy drift per-component (each |value| < 3e-5)
  5. GPU drift per-component (should not be worse than reference)
  6. Structure factor CPU/GPU agreement at FINAL stage (within 1e-4)

Usage:
    python3 score.py <test_output> <reference_output>

Exit codes:
    0 = pass (prints JSON with timing)
    1 = correctness failure
    2 = parse error
"""

import sys
import re
import json


# ===================================================================
# Thresholds
# ===================================================================
ENERGY_DRIFT_THRESHOLD = 3e-5       # Per-component energy drift max (FAIL)
ENERGY_DRIFT_IDEAL = 1e-5           # Per-component energy drift warn
GPU_DRIFT_TOLERANCE = 1e-5          # How much WORSE than reference is allowed
SF_CPU_GPU_THRESHOLD = 1e-4         # Structure factor CPU vs GPU max diff
ENERGY_MATCH_REL_TOL = 1e-5         # Relative tolerance for energy comparison
ENERGY_MATCH_ABS_TOL = 1e-3         # Absolute tolerance for energy comparison


# ===================================================================
# Parsing
# ===================================================================

def parse_output(filepath):
    """Extract all key metrics from gRASPA output."""
    with open(filepath, 'r') as f:
        text = f.read()

    result = {}

    # --- Timing ---
    m = re.search(r'Work took ([\d.]+) seconds', text)
    result['time_seconds'] = float(m.group(1)) if m else None

    # --- PseudoAtom counts at end of simulation ---
    pseudo_atoms = {}
    for m in re.finditer(r'PseudoAtom Type: (\S+)\[(\d+)\], #: (\d+)', text):
        pseudo_atoms[m.group(1)] = int(m.group(3))
    result['pseudo_atoms'] = pseudo_atoms

    # --- Component molecule counts from last cycle line ---
    last_cycle_match = None
    for m in re.finditer(
        r'(?:INITIALIZATION|EQUILIBRATION|PRODUCTION) Cycle: \d+,.+', text
    ):
        last_cycle_match = m.group(0)

    if last_cycle_match:
        m2 = re.match(r'.*?(\d+) Adsorbate Molecules', last_cycle_match)
        if m2:
            result['final_adsorbate_count'] = int(m2.group(1))
        comp_counts = {}
        for m2 in re.finditer(r'Component (\d+) \[([^\]]+)\], (\d+) Molecules', last_cycle_match):
            comp_counts[f"comp_{m2.group(1)}_{m2.group(2)}"] = int(m2.group(3))
        result['component_counts'] = comp_counts

    # --- Move statistics (last occurrence = end-of-simulation) ---
    moves = {}
    for m in re.finditer(r'Translation Performed:\s+(\d+)', text):
        moves['translation_total'] = int(m.group(1))
    for m in re.finditer(r'Translation Accepted:\s+(\d+)', text):
        moves['translation_accepted'] = int(m.group(1))
    for m in re.finditer(r'Rotation Performed:\s+(\d+)', text):
        moves['rotation_total'] = int(m.group(1))
    for m in re.finditer(r'Rotation Accepted:\s+(\d+)', text):
        moves['rotation_accepted'] = int(m.group(1))
    for m in re.finditer(r'Insertion Performed:\s+(\d+)', text):
        moves['insertion_total'] = int(m.group(1))
    for m in re.finditer(r'Insertion Accepted:\s+(\d+)', text):
        moves['insertion_accepted'] = int(m.group(1))
    for m in re.finditer(r'Deletion Performed:\s+(\d+)', text):
        moves['deletion_total'] = int(m.group(1))
    for m in re.finditer(r'Deletion Accepted:\s+(\d+)', text):
        moves['deletion_accepted'] = int(m.group(1))
    for m in re.finditer(r'Reinsertion Performed:\s+(\d+)', text):
        moves['reinsertion_total'] = int(m.group(1))
    for m in re.finditer(r'Reinsertion Accepted:\s+(\d+)', text):
        moves['reinsertion_accepted'] = int(m.group(1))
    result['moves'] = moves

    # --- Energy blocks ---
    result['final_energies'] = parse_energy_block(
        text, r'\*\*\* FINAL STAGE \*\*\*')
    result['energy_drift'] = parse_energy_block(
        text, r'\*\*\* ENERGY DRIFT \(CPU FINAL - RUNNING FINAL\) \*\*\*')
    result['gpu_drift'] = parse_energy_block(
        text, r'\*\*\* GPU DRIFT \(GPU FINAL - CPU FINAL\) \*\*\*')

    # --- Structure factor CPU/GPU agreement at FINAL stage ---
    # The SFs appear BETWEEN the two "CALCULATING FINAL STAGE ENERGY" markers
    # (the opening line and the "DONE" line both contain that string).
    # Use a targeted regex instead of split.
    final_sf_section = re.search(
        r'CALCULATING FINAL STAGE ENERGY.*?'
        r'CHECKING StructureFactors.*?'
        r'(StructureFactor.*?)(?:CHECKING Framework|VDW \+ Real on the GPU)',
        text, re.DOTALL
    )
    if final_sf_section:
        result['structure_factor_diffs'] = parse_structure_factors(final_sf_section.group(1))
    else:
        # Fallback: find the LAST block of StructureFactor lines
        all_sf_blocks = list(re.finditer(
            r'CHECKING StructureFactors.*?\n((?:StructureFactor \d+.*?\n)+)',
            text
        ))
        if all_sf_blocks:
            result['structure_factor_diffs'] = parse_structure_factors(all_sf_blocks[-1].group(1))
        else:
            result['structure_factor_diffs'] = []

    return result


def parse_energy_block(text, header_pattern):
    """Parse an energy block identified by a header pattern."""
    energies = {}
    block_match = re.search(
        header_pattern + r'.*?={5,}(.*?)={5,}', text, re.DOTALL
    )
    if not block_match:
        return energies
    block_text = block_match.group(1)
    patterns = [
        ('vdw_hh',   r"VDW \[Host-Host\]:\s+([-\d.eE+]+)"),
        ('vdw_hg',   r"VDW \[Host-Guest\]:\s+([-\d.eE+]+)"),
        ('vdw_gg',   r"VDW \[Guest-Guest\]:\s+([-\d.eE+]+)"),
        ('real_hh',  r"Real Coulomb \[Host-Host\]:\s+([-\d.eE+]+)"),
        ('real_hg',  r"Real Coulomb \[Host-Guest\]:\s+([-\d.eE+]+)"),
        ('real_gg',  r"Real Coulomb \[Guest-Guest\]:\s+([-\d.eE+]+)"),
        ('ewald_hh', r"Ewald \[Host-Host\]:\s+([-\d.eE+]+)"),
        ('ewald_hg', r"Ewald \[Host-Guest\]:\s+([-\d.eE+]+)"),
        ('ewald_gg', r"Ewald \[Guest-Guest\]:\s+([-\d.eE+]+)"),
        ('tail',     r"Tail Correction Energy:\s+([-\d.eE+]+)"),
        ('dnn',      r"DNN Energy:\s+([-\d.eE+]+)"),
        ('total',    r"Total Energy:\s+([-\d.eE+]+)"),
    ]
    for name, pattern in patterns:
        m = re.search(pattern, block_text)
        if m:
            energies[name] = float(m.group(1))
    return energies


def parse_structure_factors(text):
    """Parse StructureFactor CPU vs GPU lines. Returns list of max diffs."""
    diffs = []
    for m in re.finditer(
        r'StructureFactor \d+, CPU:\s+([-\d.]+)\s+([-\d.]+),\s+GPU:\s+([-\d.]+)\s+([-\d.]+)',
        text
    ):
        cpu_r, cpu_i = float(m.group(1)), float(m.group(2))
        gpu_r, gpu_i = float(m.group(3)), float(m.group(4))
        diffs.append(max(abs(cpu_r - gpu_r), abs(cpu_i - gpu_i)))
    return diffs


# ===================================================================
# Comparison
# ===================================================================

def values_match(test_val, ref_val):
    """Check if two numeric values match within tolerance."""
    if ref_val == 0.0 and test_val == 0.0:
        return True
    if ref_val == 0.0:
        return abs(test_val) < ENERGY_MATCH_ABS_TOL
    rel_err = abs(test_val - ref_val) / max(abs(ref_val), 1e-15)
    return rel_err < ENERGY_MATCH_REL_TOL or abs(test_val - ref_val) < ENERGY_MATCH_ABS_TOL


def compare(test, ref):
    """Compare test output against reference.
    Returns (pass, failures, warnings)."""
    failures = []
    warnings = []

    # ── CHECK 1: Move statistics — EXACT match ──
    if test.get('moves') and ref.get('moves'):
        for name, ref_val in ref['moves'].items():
            test_val = test['moves'].get(name, -1)
            if test_val != ref_val:
                failures.append(f"Move {name}: got {test_val}, expected {ref_val}")
    elif ref.get('moves'):
        failures.append("Move statistics missing from test output")

    # ── CHECK 2: PseudoAtom counts — EXACT match ──
    if test.get('pseudo_atoms') and ref.get('pseudo_atoms'):
        for atype, ref_n in ref['pseudo_atoms'].items():
            test_n = test['pseudo_atoms'].get(atype, -1)
            if test_n != ref_n:
                failures.append(f"PseudoAtom {atype}: got {test_n}, expected {ref_n}")
    elif ref.get('pseudo_atoms'):
        failures.append("PseudoAtom counts missing from test output")

    if test.get('component_counts') and ref.get('component_counts'):
        for comp, ref_n in ref['component_counts'].items():
            test_n = test['component_counts'].get(comp, -1)
            if test_n != ref_n:
                failures.append(f"Component {comp}: got {test_n}, expected {ref_n}")

    if ref.get('final_adsorbate_count') is not None:
        if test.get('final_adsorbate_count') != ref['final_adsorbate_count']:
            failures.append(
                f"Final adsorbate count: got {test.get('final_adsorbate_count')}, "
                f"expected {ref['final_adsorbate_count']}")

    # ── CHECK 3: Final energies — within tolerance ──
    if test.get('final_energies') and ref.get('final_energies'):
        for name, ref_val in ref['final_energies'].items():
            test_val = test['final_energies'].get(name)
            if test_val is None:
                failures.append(f"Energy {name}: missing in test output")
            elif not values_match(test_val, ref_val):
                failures.append(
                    f"Energy {name}: got {test_val}, expected {ref_val}, "
                    f"diff={abs(test_val - ref_val):.6e}")

    # ── CHECK 4: Energy drift — each component |val| < 3e-5 ──
    if test.get('energy_drift'):
        for name, val in test['energy_drift'].items():
            av = abs(val)
            if av > ENERGY_DRIFT_THRESHOLD:
                failures.append(
                    f"Energy drift [{name}]: {val:.6e} exceeds {ENERGY_DRIFT_THRESHOLD}")
            elif av > ENERGY_DRIFT_IDEAL:
                warnings.append(
                    f"Energy drift [{name}]: {val:.6e} above ideal {ENERGY_DRIFT_IDEAL}")
    else:
        warnings.append("Energy drift block not found in test output")

    # ── CHECK 5: GPU drift — should not be WORSE than reference ──
    # GPU drift is inherent CPU/GPU numerical difference. We compare against
    # the reference's GPU drift: if the test is no worse, it's fine.
    if test.get('gpu_drift') and ref.get('gpu_drift'):
        for name, test_val in test['gpu_drift'].items():
            ref_val = ref['gpu_drift'].get(name, 0.0)
            # Fail only if test drift is significantly worse than reference
            if abs(test_val) > abs(ref_val) + GPU_DRIFT_TOLERANCE:
                failures.append(
                    f"GPU drift [{name}]: {test_val:.6e} worse than reference {ref_val:.6e} "
                    f"(tolerance {GPU_DRIFT_TOLERANCE})")
    elif test.get('gpu_drift'):
        # No reference GPU drift to compare against — use absolute threshold
        for name, val in test['gpu_drift'].items():
            if abs(val) > ENERGY_DRIFT_THRESHOLD:
                failures.append(
                    f"GPU drift [{name}]: {val:.6e} exceeds {ENERGY_DRIFT_THRESHOLD}")

    # ── CHECK 6: Structure factor CPU/GPU agreement ──
    if test.get('structure_factor_diffs'):
        max_diff = max(test['structure_factor_diffs'])
        n = len(test['structure_factor_diffs'])
        if max_diff > SF_CPU_GPU_THRESHOLD:
            failures.append(
                f"Structure factor CPU/GPU mismatch: max diff = {max_diff:.6e} "
                f"(across {n} SFs), threshold = {SF_CPU_GPU_THRESHOLD}")
    else:
        warnings.append("Structure factor CPU/GPU check not found in test output")

    return len(failures) == 0, failures, warnings


# ===================================================================
# Main
# ===================================================================

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <test_output> <reference_output>", file=sys.stderr)
        sys.exit(2)

    try:
        test = parse_output(sys.argv[1])
        ref = parse_output(sys.argv[2])
    except Exception as e:
        print(json.dumps({"status": "parse_error", "error": str(e)}))
        sys.exit(2)

    if test['time_seconds'] is None:
        print(json.dumps({"status": "parse_error",
                          "error": "Could not find 'Work took X seconds' in test output"}))
        sys.exit(2)

    passed, failures, warnings = compare(test, ref)

    # Summaries
    drift_summary = {}
    for key in ('energy_drift', 'gpu_drift'):
        if test.get(key):
            nonzero = {k: f"{v:.6e}" for k, v in test[key].items() if abs(v) > 1e-10}
            if nonzero:
                drift_summary[key] = nonzero

    sf_max = max(test['structure_factor_diffs']) if test.get('structure_factor_diffs') else None

    result = {
        "status": "pass" if passed else "fail",
        "time_seconds": test['time_seconds'],
        "baseline_seconds": ref['time_seconds'],
        "speedup": (round(ref['time_seconds'] / test['time_seconds'], 4)
                    if test['time_seconds'] and ref['time_seconds'] and test['time_seconds'] > 0
                    else None),
        "failures": failures if not passed else [],
        "warnings": warnings,
        "final_adsorbate_count": test.get('final_adsorbate_count'),
        "total_energy": test.get('final_energies', {}).get('total'),
        "drift_summary": drift_summary,
        "sf_max_cpu_gpu_diff": sf_max,
    }

    print(json.dumps(result, indent=2))
    sys.exit(0 if passed else 1)


if __name__ == '__main__':
    main()
