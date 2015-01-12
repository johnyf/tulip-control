"""Tests for transys.transys (part of transys subpackage)"""
import logging
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)
logging.getLogger('tulip.transys.products').setLevel(logging.DEBUG)
from tulip import transys as trs
from tulip.transys import products


def ts_test():
    ts = trs.TransitionSystem()
    ts.add_nodes_from([0, 1, 2])
    ts.initial_nodes.update([0, 2])
    ts.vars = {'a': 'boolean', 'b': (0, 3)}
    ts.env_vars = {'a'}
    ts.add_node(0, a=True, b=1)
    ts.add_node(1, a=False, b=3)
    ts.add_node(2, a=True, b=2)
    ts.add_edge(0, 1, **{"a'": False})
    ts.add_edge(0, 1, **{"a'": True, "b'": 2})
    assert ts.is_consistent()


def ba_ts_prod_test():
    ts = ts_test()
    ba = ba_test()

    ba_ts = products.ba_ts_sync_prod(ba, ts)
    check_prodba(ba_ts)

    (ts_ba, persistent) = products.ts_ba_sync_prod(ts, ba)

    states = {('s0', 'q1'), ('s1', 'q0'),
              ('s2', 'q0'), ('s3', 'q0')}

    assert set(ts_ba.states) == states
    assert persistent == {('s0', 'q1')}

    ba_ts.save('prod.pdf')
    return ba_ts


def check_prodba(ba_ts):
    states = {('s0', 'q1'), ('s1', 'q0'),
              ('s2', 'q0'), ('s3', 'q0')}

    assert set(ba_ts.states) == states

    assert set(ba_ts.states.initial) == {('s0', 'q1'), ('s1', 'q0')}

    assert (
        ba_ts.transitions.find(
            [('s0', 'q1')], [('s1', 'q0')]
        )[0][2]['guard'] == set())

    for si, sj in [('s1', 's2'), ('s2', 's3')]:
        assert (
            ba_ts.transitions.find(
                [(si, 'q0')], [(sj, 'q0')]
            )[0][2]['guard'] == set())

    assert (
        ba_ts.transitions.find(
            [('s3', 'q0')], [('s0', 'q1')]
        )[0][2]['guard'] == {'p'})
