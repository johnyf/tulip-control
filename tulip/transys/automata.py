# Copyright (c) 2013-2014 by California Institute of Technology
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the California Institute of Technology nor
#    the names of its contributors may be used to endorse or promote
#    products derived from this software without specific prior
#    written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CALTECH
# OR THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
"""Automata Module"""
from __future__ import absolute_import
import logging
logger = logging.getLogger(__name__)
import copy
from collections import Iterable
from pprint import pformat
from tulip.transys.labeled_graphs import (
    LabeledDiGraph, str2singleton, prepend_with)
from tulip.transys.mathset import SubSet, PowerSet
from tulip.transys.transys import GameGraph





def tuple2ba(S, S0, Sa, Sigma_or_AP, trans, name='ba', prepend_str=None,
             atomic_proposition_based=True):
    """Create a Buchi Automaton from a tuple of fields.

    defines Buchi Automaton by a tuple (S, S0, Sa, \\Sigma, trans)
    (maybe replacing \\Sigma by AP since it is an AP-based BA ?)

    See Also
    ========
    L{tuple2fts}

    @param S: set of states
    @param S0: set of initial states, must be \\subset S
    @param Sa: set of accepting states
    @param Sigma_or_AP: Sigma = alphabet
    @param trans: transition relation, represented by list of triples::
            [(from_state, to_state, guard), ...]
    where guard \\in \\Sigma.

    @param name: used for file export
    @type name: str

    @rtype: L{BuchiAutomaton}
    """
    # args
    if not isinstance(S, Iterable):
        raise TypeError('States S must be iterable, even for single state.')
    if not isinstance(S0, Iterable) or isinstance(S0, str):
        S0 = [S0]
    if not isinstance(Sa, Iterable) or isinstance(Sa, str):
        Sa = [Sa]
    # comprehensive names
    states = S
    initial_states = S0
    accepting_states = Sa
    alphabet_or_ap = Sigma_or_AP
    transitions = trans
    # prepending states with given str
    if prepend_str:
        logger.debug('Given string:\n\t' + str(prepend_str) + '\n' +
                     'will be prepended to all states.')
    states = prepend_with(states, prepend_str)
    initial_states = prepend_with(initial_states, prepend_str)
    accepting_states = prepend_with(accepting_states, prepend_str)

    ba = BuchiAutomaton(atomic_proposition_based=atomic_proposition_based)
    ba.name = name

    ba.states.add_from(states)
    ba.states.initial |= initial_states
    ba.states.accepting |= accepting_states

    if atomic_proposition_based:
        ba.alphabet.math_set |= alphabet_or_ap
    else:
        ba.alphabet.add(alphabet_or_ap)
    for transition in transitions:
        (from_state, to_state, guard) = transition
        [from_state, to_state] = prepend_with([from_state, to_state],
                                              prepend_str)
        # convention
        if atomic_proposition_based:
            if guard is None:
                guard = set()
            guard = str2singleton(guard)
        ba.transitions.add(from_state, to_state, letter=guard)
    return ba




class ParityGame(GameGraph):
    """GameGraph equipped with coloring.

    Define as C{k} the highest color that
    occurs infinitely many times.

    If C{k} is even, then Player 0 wins.
    Otherwise Player 1 wins (C{k} is odd).
    So the winner is Player (k mod 2).

    To define the number of colors C{c}:

    >>> p = ParityGame(c=4)

    Note that the colors are: 0, 1, ..., c-1

    See also
    ========
    L{transys.GameGraph}
    """

    def __init__(self, c=2):
        node_label_types = [{
            'name': 'color',
            'values': range(c),
            'default': 0}]
        super(ParityGame, self).__init__(node_label_types, [])

    def __str__(self):
        s = (
            'Parity Game\n'
            '-----------\n'
            'n: node, p: player, c: color\n\n')
        for node, attr in self.states(data=True):
            s += 'nd = {node}, p = {player}, c = {color}\n'.format(
                npde=node, player=attr['player'], color=attr['color'])
        s += '\n{t}'.format(t=self.transitions)
        return s

    @property
    def max_color(self):
        max_c = -1
        # node = None
        for x in self:
            if self.node[x]['color'] > max_c:
                max_c = self.node[x]['color']
                # node = x
        return max_c
