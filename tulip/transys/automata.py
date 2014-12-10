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
from collections import Iterable
from pprint import pformat
from tulip.transys.labeled_graphs import (
    LabeledDiGraph, str2singleton, prepend_with)
from tulip.transys.transys import GameGraph


class Automaton(LabeledDiGraph):
    """Alternating acceptor of (in)finite trees.

    An acceptor represents an indicator function of a set.
    The set may contain trees (unary trees are words).

    From all the possible paths in the acceptor's graph,
    some correspond to elements of the set it describes,
    but others do not belong to this set.

    In order to distinguish elements that can be generated by
    the graph, but are not contained in the set,
    two attributes are used:

      - path quantification
      - acceptance condition


    Attributes
    ==========
    Besides initial nodes
    (from [`LabeledDiGraph`] - to be split and renamed),
    an [`Automaton`] has the attributes:

      - `"universal_nodes"` is a subset of nodes.
        Remaining nodes are existentially quantified.


      - `"acceptance"` that can be:

        - `"finite"`, signifying finite words or trees
        - `"Buchi"` (or `"weak Buchi"`)
        - `"coBuchi"`
        - `"generalized Buchi"`
        - `"Rabin"`
        - `"Streett"`
        - `"Muller"`
        - `"parity"`


      - `"accepting_sets"` that contains the sets of nodes defining
        the acceptance condition:

        - `set` of nodes for (weak/co-) Buchi and finite
          For co-Buchi this is the avoidance set (FG!)
        - `list` of `dict`s of `set`s of nodes for Rabin and Streett
          The keys should be: `"GF"` and `"FG!"`,
        - `list` of `set`s of nodes for Muller and generalized Buchi,
        - a 2-`tuple` of (min, max) colors for parity
          The coloring function is defined by directly annotating the nodes.

        For convenience, the above are initialized so that
        the automaton represent an empty language.
        The user is responsible for adhering to the above conventions.

        To check compliance to the above, call [`Automaton.check_sanity`].

        Deprecated: Note that by using `SubSet(automaton)` for each set,
        you can ensure consistency if more nodes are added later.
        Alternatively, write a function that checks consistency of
        each acceptance condition wrt the automaton's nodes.


      - `alphabet`: `dict` mapping from symbols to domains.
        Invariants of existential nodes are models or
        formulae over `alphabet`.


      - `directions`: Like `alphabet`, but for universal nodes.
        If `directions` is `None`, then the automaton recognizes words.


      - `guards` defines the representation of edge labels and can be:

        - `"formula"` meaning that each edge is annotated with
          a Boolean formula as `str` or AST. (Typically a conjunction.)
        - `"enumeration"` meaning that each edge is annotated with a letter,
          as `set` or `dict` (closed world semantics).


    Related concepts
    ================
    Automata represent languages in a way suitable for
    testing if a given trace is a member of the language.
    The represented language is not readily accessible,
    because its generation requires solving a search problem.

    If you want to represent a language constructively,
    then use a [`KripkeStructure`] instead.
    That is equivalent to a universal Buchi acceptor
    whose nodes are all accepting.

    To represent transductions use a [`Transducer`].


    Remark
    ======
    In `ltl2dstar` documentation L denotes a "good" set.
    To avoid ambiguity, `dict`s with explicit modalities are used as labels.


    References
    ==========
    - Def. 10.53, p.801, U{[BK08]
      <http://tulip-control.sourceforge.net/doc/bibliography.html#bk08>}
    - U{ltl2dstar<http://ltl2dstar.de/>} documentation


    See also
    ========
    [`KripkeStructure`], [`Transducer`]
    """

    def __init__(self, acceptance='Buchi', alphabet=None,
                 directions=None, universal_nodes=None, guards='boolean'):
        if universal_nodes is None:
            universal_nodes = set()
        if alphabet is None:
            alphabet = dict()
        # init attributes
        self.universal_nodes = universal_nodes
        self.acceptance = acceptance
        self.accepting_sets = self._init_accepting_sets(acceptance)
        # default: powerset alphabet
        alphabet = alphabet
        self.alphabet = alphabet
        self.directions = directions
        self.guards = guards
        # type checking for edge labeling
        edge_label_types = [
            {'name': 'guard',
             'values': None,
             'setter': True}]
        super(Automaton, self).__init__(edge_label_types=edge_label_types)
        # formatting options
        self._transition_dot_label_format = {
            'guard': '',
            'type?label': '',
            'separator': '\n'}
        self._transition_dot_mask = dict()
        self.dot_node_shape = {'normal': 'circle',
                               'accepting': 'doublecircle'}
        self.default_export_fname = 'fsa'

    def __str__(self):
        show_node_data = (self.acceptance == 'parity')
        word_tree = (self.directions is None or len(self.directions) <= 1)
        f = lambda x: pformat(x, indent=3)
        s = (
            '{hl}\n Alternating {acceptance} {word_tree} automaton\n'
            '{hl}\n'
            'Alphabet:\n'
            '{self.alphabet}\n\n'
            'Directions:\n'
            '{self.directions}\n\n'
            'Nodes:\n'
            '{nodes}\n\n'
            'Initial nodes:\n'
            '{init_nodes}\n\n'
            'Universal nodes:\n'
            '{universal_nodes}\n\n'
            'Existential nodes: the rest\n\n'
            'Accepting sets:\n'
            '{accepting_sets}\n\n'
            'Edges with guards:\n'
            '{edges}\n{hl}\n').format(
                hl=40 * '-',
                acceptance=self.acceptance,
                word_tree=word_tree,
                self=self,
                nodes=f(self.nodes(data=show_node_data)),
                init_nodes=f(self.states.initial),
                universal_nodes=f(self.universal_nodes),
                accepting_sets=f(self.accepting_sets),
                edges=f(self.transitions(data=True)))
        return s

    def _init_accepting_sets(self, acceptance):
        if acceptance in {'finite', 'Buchi', 'coBuchi', 'weak Buchi'}:
            a = set()
        elif acceptance in {'Muller', 'generalized Buchi', 'Rabin', 'Streett'}:
            a = list()
        elif acceptance == 'parity':
            a = (None, None)
        else:
            raise ValueError('unknown acceptance: {s}'.format(s=acceptance))
        return a

    def check_sanity(self):
        a = self.acceptance
        s = self.accepting_sets
        f = lambda x: all(u in self for u in x)
        if a in {'finite', 'Buchi', 'coBuchi', 'weak Buchi'}:
            assert f(s)
        elif a in {'Muller', 'generalized Buchi'}:
            for x in s:
                assert f(s)
        elif a in {'Rabin', 'Streett'}:
            for x in s:
                assert len(x) == 2
                assert set(x) == {'[]<>', '<>[]!'}
                assert f(x['[]<>'])
                assert f(x['<>[]!'])
        elif a == 'parity':
            assert len(s) == 2
        else:
            raise Exception('Unknown acceptance: {a}'.format(a=a))
        return True


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
    raise NotImplementedError('currently defunct')
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

    ba = Automaton(atomic_proposition_based=atomic_proposition_based)
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
