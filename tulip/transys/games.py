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
    """A labeled game graph with parity winning condition.

    When adding states, you have to say
    which player controls the outgoing transitions.
    Use C{networkx} state labels for that:

      >>> g = GameGraph()
      >>> g.add_node('s0', player=0, color=2)


    Reference
    =========
    1. Chatterjee K.; Henzinger T.A.; Jobstmann B.
       Environment Assumptions for Synthesis
       CONCUR'08, LNCS 5201, pp. 147-161, 2008

    2. Chatterjee K.; Henzinger T.A.; Piterman N.
       Strategy Logic
       UCB/EECS-2007-78
    """

    def __init__(self):
        super(GameGraph, self).__init__()
        self.colors = None

    def __str__(self):
        s = (
            'Parity Game\n'
            '-----------\n'
            'n: node, p: player, c: color\n\n')
        for node, attr in self.nodes_iter(data=True):
            s += 'nd = {node}, p = {player}, c = {color}\n'.format(
                npde=node, player=attr['player'], color=attr['color'])
        s += '\n{t}'.format(t=self.edges())
        return s

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
        u = e[0]
        return self.node[u]['player']

    def to_pydot(self):
        pass

    @property
    def max_color(self):
        c = -1
        for u, d in self.nodes_iter(data=True):
            c = max(c, d['color'])
        return c
