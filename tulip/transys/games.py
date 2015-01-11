# Copyright (c) 2014 by California Institute of Technology
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
"""Two-player infinite games on finite graphs."""
import logging
logger = logging.getLogger(__name__)
import networkx as nx


class GameGraph(nx.MultiDiGraph):
    """Store a game graph.

    When adding states, you have to say
    which player controls the outgoing transitions.
    Use C{networkx} state labels for that:

      >>> g = GameGraph()
      >>> g.states.add('s0', player=0)

    See also
    ========
    L{automata.ParityGame}

    Reference
    =========
    1. Chatterjee K.; Henzinger T.A.; Jobstmann B.
       Environment Assumptions for Synthesis
       CONCUR'08, LNCS 5201, pp. 147-161, 2008
    """

    def __init__(self, node_label_types, edge_label_types):
        node_label_types += [{
            'name': 'player',
            'values': {0, 1},
            'default': 0}]
        super(GameGraph, self).__init__(node_label_types,
                                        edge_label_types)

    def player_states(self, n):
        """Return states controlled by player C{n}.

        'controlled' means that player C{n}
        gets to decide the successor state.

        @param n: player index (id number)
        @type n: 0 or 1

        @return: set of states
        @rtype: C{set}
        """
        return {x for x in self if self.node[x]['player'] == n}

    def edge_controlled_by(self, e):
        """Return the index of the player controlling edge C{e}.

        @type e: 2-tuple of nodes C{(n1, n2)}

        @rtype: integer 0 or 1
        """
        from_state = e[0]
        return self.node[from_state]['player']


class LabeledGameGraph(GameGraph):
    """Game graph with labeled states.

    Its contraction is a Kripke structure.
    Given a Kripke structure and a partition of propositions,
    then the corresponding labeled game graph
    can be obtained by graph expansion.

    Reference
    =========
    1. Chatterjee K.; Henzinger T.A.; Piterman N.
       Strategy Logic
       UCB/EECS-2007-78
    """

    def __init__(self):
        ap_labels = PowerSet()
        node_label_types = [
            {'name': 'ap',
             'values': ap_labels,
             'setter': ap_labels.math_set,
             'default': set()}]
        super(LabeledGameGraph, self).__init__(node_label_types)
        self.atomic_propositions = self.ap
        # dot formatting
        self._state_dot_label_format = {
            'ap': '',
            'type?label': '',
            'separator': '\n'}
        self.dot_node_shape = {'normal': 'rectangle'}


class ParityGame(GameGraph):
    """GameGraph with coloring.

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
