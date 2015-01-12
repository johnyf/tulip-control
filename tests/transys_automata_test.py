#!/usr/bin/env python
"""Tests for transys.automata (part of transys subpackage)"""
from tulip import transys as trs


def rabin_test():
    dra = trs.Automaton(acceptance='Rabin')
    print(dra)

    dra.add_nodes_from(xrange(10))
    dra.accepting_sets.append({'[]<>': {1, 2}, '<>[]!': {3, 4}})

    assert isinstance(dra.accepting_sets, list)
    assert dra.accepting_sets[0]['[]<>'] == {1, 2}
    assert dra.accepting_sets[0]['<>[]!'] == {3, 4}
    assert dra.check_sanity()

    dra = trs.Automaton(acceptance='Rabin')
    dra.add_nodes_from(xrange(5))
    dra.accepting_sets.append({'[]<>': [1], '<>[]!': [2, 4]})
    dra.accepting_sets[0]['[]<>'].append(2)
    dra.accepting_sets.append({'[]<>': [2, 1], '<>[]!': []})

    assert isinstance(dra.accepting_sets, list)
    s = dra.accepting_sets[0]
    assert set(s['[]<>']) == {1, 2}
    assert set(s['<>[]!']) == {2, 4}
    s = dra.accepting_sets[1]
    assert set(s['[]<>']) == {2, 1}
    assert set(s['<>[]!']) == set()
    assert dra.check_sanity()
