"""Checks the IdentitySwap-IdealGas example against its exact answer.

The example is a two-component ideal gas, so <N_c> = beta * f_c * V exactly
(see IdentitySwap-IdealGas/README.md). The loading test is skipped until the
example has been run and has an output.txt; the parser itself is exercised on
the shipped XeKr-Mixture output so that it stays in step with the output format.
"""
import os
import re

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
FOLDER = os.path.join(HERE, "IdentitySwap-IdealGas")
OUTPUT = os.path.join(FOLDER, "output.txt")

# beta * f_c * V for the two components of IdentitySwap-IdealGas/simulation.input
EXACT_LOADING = {"A": 0.3, "B": 2.7}
# |<N_c> - B_c| must not exceed max(ErrorBar_c, REL_FLOOR * B_c). The old
# pair-selection rule misses by about 30% (A) and 12% (B).
REL_FLOOR = 0.03

LOADING_HEADER = "BLOCK AVERAGES (LOADING: # MOLECULES)"
COMPONENT_LINE = re.compile(r"COMPONENT \[(\d+)\] \((.+)\)")
OVERALL_LINE = re.compile(r"Overall: Average:\s*([-+0-9.eE]+),\s*ErrorBar:\s*([-+0-9.eE]+)")


def overall_loadings(text):
    """{component name: (average, error bar)} from the last loading block of a gRASPA output.

    Each component prints two 'Overall' lines (absolute loading, then excess
    loading); the first one is taken."""
    start = text.rfind(LOADING_HEADER)
    if start < 0:
        raise ValueError("no '%s' block in the output" % LOADING_HEADER)
    result = {}
    component = None
    for line in text[start + len(LOADING_HEADER):].splitlines():
        if "BLOCK AVERAGES" in line:   # the next kind of block average
            break
        m = COMPONENT_LINE.match(line.strip())
        if m:
            component = m.group(2)
            continue
        m = OVERALL_LINE.match(line.strip())
        if m and component is not None and component not in result:
            result[component] = (float(m.group(1)), float(m.group(2)))
    return result


def test_loading_parser_on_shipped_xekr_output():
    with open(os.path.join(HERE, "XeKr-Mixture", "output.txt")) as f:
        loadings = overall_loadings(f.read())
    assert set(loadings) >= {"Kr", "Xe"}, loadings
    assert loadings["Xe"][0] > loadings["Kr"][0] > 0, loadings


@pytest.mark.skipif(not os.path.exists(OUTPUT),
                    reason="IdentitySwap-IdealGas has not been run (no output.txt)")
def test_identity_swap_ideal_gas_loadings():
    with open(OUTPUT) as f:
        loadings = overall_loadings(f.read())
    for name, exact in EXACT_LOADING.items():
        assert name in loadings, "component %s not found in %s: %s" % (name, OUTPUT, sorted(loadings))
        mean, errorbar = loadings[name]
        tolerance = max(errorbar, REL_FLOOR * exact)
        print("IdentitySwap-IdealGas %s: <N> = %.5f +- %.5f, exact %.5f, tolerance %.5f"
              % (name, mean, errorbar, exact, tolerance))
        assert abs(mean - exact) <= tolerance, (
            "<N_%s> = %.5f differs from the exact %.5f by %.5f (tolerance %.5f); "
            "the identity-swap proposal or acceptance is biased" % (name, mean, exact, abs(mean - exact), tolerance))
