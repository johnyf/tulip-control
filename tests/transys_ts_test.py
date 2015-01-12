"""Tests for transys.transys (part of transys subpackage)"""
import logging
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)
logging.getLogger('tulip.transys.products').setLevel(logging.DEBUG)
from tulip import transys as trs
from tulip.transys import products


def ts_test():
    ts = trs.TransitionSystem()

    ts.states.add('s0')
    assert 's0' in ts
    assert 's0' in ts.node
    assert 's0' in ts.states

    states = {'s0', 's1', 's2', 's3'}
    ts.states.add_from(states)
    assert set(ts.states) == states

    ts.transitions.add('s0', 's1')
    ts.transitions.add_from([('s1', 's2'), ('s2', 's3'), ('s3', 's0')])

    ts.states.initial.add('s0')
    ts.states.initial.add_from({'s0', 's1'})

    ts.atomic_propositions.add('p')
    assert set(ts.atomic_propositions) == {'p'}

    ts.states['s0']['ap'] = {'p'}
    ts.states['s1']['ap'] = set()
    ts.states['s2']['ap'] = set()
    ts.states['s3']['ap'] = set()

    assert ts.states['s0']['ap'] == {'p'}
    for state in {'s1', 's2', 's3'}:
        assert ts.states[state]['ap'] == set()

    logger.debug(ts)
    return ts




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
