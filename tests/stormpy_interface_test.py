# This test includes integration tests for the interfaces.stormpy module,
# using the example from the EECI 2020 computer lab.
# Details of these tests can be found at
# the Stormpy computer session of EECI 2020:
#
# https://www.cds.caltech.edu/~murray/wiki/index.php?title=EECI_2020:_Computer_Session:_Stormpy
#
# All the tests here are run only if stormpy is available.
#
# * light_analysis_test() corresponds to Example 1: traffic light (analysis)
# * compose_test() corresponds to Example 2: traffic light + vehicle (analysis)
# * synthesis_test() corresponds to Example 4: policy synthesis

import os
import unittest
from tulip.transys.compositions import synchronous_parallel, apply_policy

# Set up paths to all the models
model_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "models")
light_path = os.path.join(model_path, "light.pm")
ma_path = os.path.join(model_path, "ma.nm")
mh_path = os.path.join(model_path, "mh.pm")
out_model_path = os.path.join(model_path, "tmpout.nm")

skip_test = False
try:
    # import stormpy here so that the import error is catched in python2.7
    import stormpy
    from tulip.interfaces import stormpy as stormpy_int
except ImportError:
    skip_test = True


@unittest.skipIf(skip_test, "stormpy not available")
def light_analysis_test():
    # Build models from prism files
    light = stormpy_int.to_tulip_transys(light_path)

    # Check properties
    formula = 'P=? [ F ("green") ]'

    # Model checking
    result = stormpy_int.model_checking(light, formula, out_model_path)
    os.remove(out_model_path)

    assert len(light.states) == 2
    expected_result = [
        {"labels": {"green"}, "result": 1.0},
        {"labels": {"red"}, "result": 1.0},
    ]

    _check_result(light, result, expected_result)


@unittest.skipIf(skip_test, "stormpy not available")
def compose_test():
    # Build models from prism files
    mh = stormpy_int.to_tulip_transys(mh_path)
    light = stormpy_int.to_tulip_transys(light_path)

    # Compose models
    composed = synchronous_parallel([mh, light])

    # Check properties
    formula = 'P=? [ "green" U "h6" ]'

    # Model checking
    result = stormpy_int.model_checking(composed, formula, out_model_path)
    os.remove(out_model_path)

    assert len(composed.states) == 14
    expected_result = [
        {"labels": {"h6", "green"}, "result": 1.0},
        {"labels": {"red", "h6"}, "result": 1.0},
        {"labels": {"h2", "green"}, "result": 0.42122366709344317},
        {"labels": {"red", "h2"}, "result": 0.0},
        {"labels": {"h4", "green"}, "result": 0.7256235827664395},
        {"labels": {"red", "h4"}, "result": 0.0},
        {"labels": {"h5", "green"}, "result": 0.9523809523809521},
        {"labels": {"red", "h5"}, "result": 0.0},
        {"labels": {"green", "h1"}, "result": 0.3209323177854804},
        {"labels": {"red", "h1"}, "result": 0.0},
        {"labels": {"h0", "green"}, "result": 0.2445198611698898},
        {"labels": {"red", "h0"}, "result": 0.0},
        {"labels": {"green", "h3"}, "result": 0.5528560630601442},
        {"labels": {"red", "h3"}, "result": 0.0},
    ]

    _check_result(composed, result, expected_result)


@unittest.skipIf(skip_test, "stormpy not available")
def synthesis_test():
    # Build models from prism files
    ma = stormpy_int.to_tulip_transys(ma_path)
    mh = stormpy_int.to_tulip_transys(mh_path)
    light = stormpy_int.to_tulip_transys(light_path)

    # Compose models
    composed = synchronous_parallel([ma, mh, light])

    # Check properties
    safety = '!("h4" & "a4") & !("red" & ("a8" | "a4"))'
    reach = '"a9"'
    formula = "Pmax=? [ ({}) U ({}) ]".format(safety, reach)

    # Construct policy
    (result, policy) = stormpy_int.model_checking(
        composed, formula, out_model_path, True
    )
    os.remove(out_model_path)

    assert abs(result[list(composed.states.initial)[0]] - 0.7429340826573935) < 1e-6

    # Get the MC induced by applying policy_opt on model
    induced_mc = apply_policy(composed, policy)

    # Model checking
    result = stormpy_int.model_checking(induced_mc, formula, out_model_path)
    os.remove(out_model_path)

    expected_result = [
        {"labels": {"green", "h5", "a8"}, "result": 0.745341614906832},
        {"labels": {"h5", "red", "a8"}, "result": 0.0},
        {"labels": {"green", "h0", "a8"}, "result": 0.7429340826573935},
        {"labels": {"h0", "red", "a8"}, "result": 0.0},
        {"labels": {"green", "h1", "a8"}, "result": 0.7287988161717747},
        {"labels": {"h1", "a8", "red"}, "result": 0.0},
        {"labels": {"green", "a8", "h4"}, "result": 0.6059687926071806},
        {"labels": {"red", "a8", "h4"}, "result": 0.0},
        {"labels": {"green", "a8", "h6"}, "result": 0.745341614906832},
        {"labels": {"red", "a8", "h6"}, "result": 0.0},
        {"labels": {"green", "a8", "h3"}, "result": 0.3861264940442798},
        {"labels": {"red", "a8", "h3"}, "result": 0.0},
        {"labels": {"green", "h2", "a8"}, "result": 0.6458232194557872},
        {"labels": {"h2", "a8", "red"}, "result": 0.0},
        {"labels": {"green", "h5", "a4"}, "result": 0.9523809523809523},
        {"labels": {"h5", "red", "a4"}, "result": 0.0},
        {"labels": {"green", "h0", "a4"}, "result": 0.9520897806888625},
        {"labels": {"h0", "red", "a4"}, "result": 0.0},
        {"labels": {"green", "h1", "a4"}, "result": 0.9501789664595234},
        {"labels": {"h1", "a4", "red"}, "result": 0.0},
        {"labels": {"green", "a4", "h4"}, "result": 0.0},
        {"labels": {"red", "a4", "h4"}, "result": 0.0},
        {"labels": {"green", "a4", "h6"}, "result": 0.9523809523809521},
        {"labels": {"red", "a4", "h6"}, "result": 0.0},
        {"labels": {"green", "a4", "h3"}, "result": 0.8264462809917354},
        {"labels": {"red", "a4", "h3"}, "result": 0.0},
        {"labels": {"green", "h2", "a4"}, "result": 0.9357284338501467},
        {"labels": {"h2", "a4", "red"}, "result": 0.0},
        {"labels": {"green", "h5", "a9"}, "result": 1.0},
        {"labels": {"h5", "a9", "red"}, "result": 1.0},
        {"labels": {"green", "h0", "a9"}, "result": 1.0},
        {"labels": {"h0", "a9", "red"}, "result": 1.0},
        {"labels": {"green", "a9", "h1"}, "result": 1.0},
        {"labels": {"a9", "red", "h1"}, "result": 1.0},
        {"labels": {"green", "a9", "h4"}, "result": 1.0},
        {"labels": {"a9", "red", "h4"}, "result": 1.0},
        {"labels": {"green", "a9", "h6"}, "result": 1.0},
        {"labels": {"a9", "h6", "red"}, "result": 1.0},
        {"labels": {"green", "a9", "h3"}, "result": 1.0},
        {"labels": {"a9", "h3", "red"}, "result": 1.0},
        {"labels": {"green", "a9", "h2"}, "result": 1.0},
        {"labels": {"a9", "h2", "red"}, "result": 1.0},
    ]

    _check_result(composed, result, expected_result)


def _check_result(model, result, expected_result):
    ap_set = [res["labels"] for res in expected_result]

    # Examine result
    for state in model.states:
        assert model.states[state]["ap"] in ap_set
        # Make sure that each atomic proposition belongs to only one state
        ap_set.remove(model.states[state]["ap"])
        expected_probability = [
            res["result"]
            for res in expected_result
            if res["labels"] == model.states[state]["ap"]
        ]
        assert abs(result[state] - expected_probability[0]) < 1e-6
