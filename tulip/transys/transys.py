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
"""Transition System Module"""
import logging
logger = logging.getLogger(__name__)
from pprint import pformat
import networkx as nx
from tulip.transys.labeled_graphs import SystemGraph, check_value


class TransitionSystem(SystemGraph):
    """Enumerated or semi-symbolic transition system.

    Next node
    =========
    The attribute `owner` defines who selects the next node.
    It can be assigned the values `"env"` or `"sys"`.


    Node and edge labels
    ====================
    The `vars` attribute is a `dict` mapping variable symbols to domains.
    The `env_vars` attribute is a `set` of keys from `vars`
    that are controlled by the environment.

    Node and edge attributes can be:
      - assignments
      - a formula (key: `"formula"`)

    These are over:
      - unprimed `vars` for nodes
      - primed and unprimed `vars` for edges

    An assignment is a key from `vars` (possibly primed)
    paired with a value.
    An assignment is interpreted as the corresponding conjunction.
    So a partial assignment is understood using open world semantics.
    A formula can be a `str` or syntax tree.


    Open vs closed
    ==============
    If `env_vars` is empty, then the system is closed.
    Otherwise the system is open.
    The quantification of `vars` is the
    main difference from a Kripke structure.


    See Also
    ========
    L{KripkeStructure}
    """

    def __init__(self):
        super(TransitionSystem, self).__init__()
        self.owner = 'sys'
        self.vars = dict()
        self.env_vars = set()

    def __str__(self):
        return (
            '{hl}\n'
            'Transition system:\n'
            '{hl}\n'
            'Owner:\n'
            '{self.owner}\n'
            'Variables:\n'
            '{self.vars}\n'
            'Environment variables:\n'
            '{self.env_vars}\n'
            'Nodes with assignments to vars:\n'
            '{nodes}\n'
            'Initial nodes:\n'
            '{self.initial_nodes}\n'
            'Edges with actions:\n'
            '{edges}\n'
            '{hl}\n').format(
                hl=40 * '-',
                self=self,
                nodes=_dumps_nodes(self),
                edges=pformat(self.edges(data=True), indent=3))

    def to_pydot(self):
        g = nx.MultiDiGraph()
        for u, d in self.nodes_iter(data=True):
            label = ', '.join(
                '{k} = {v}'.format(k=k, v=v)
                for k, v in d.iteritems()
                if k in self.vars)
            g.add_node(u, label=label, shape='box')
        for u, v, d in self.edges_iter(data=True):
            label = 'sys: {sys}\n env: {env}'.format(
                env=self.env_actions,
                sys=self.sys_actions)
        return nx.to_pydot(g)

    def is_consistent(self):
        """Return `True` if attributes are conformant."""
        # TODO: check type consistency of formulae
        super(TransitionSystem, self).is_consistent()
        assert self.owner in {'env', 'sys'}
        assert set(self.env_vars).issubset(self.vars)
        for u, d in self.nodes_iter(data=True):
            for k, v in d.iteritems():
                if k in self.vars:
                    check_value(v, self.vars[k])
        for u, v, d in self.edges_iter(data=True):
            for k, v in d.iteritems():
                # primed ?
                if k.endswith("'"):
                    var = k[:-1]
                else:
                    var = k
                if var in self.vars:
                    check_value(v, self.vars[var])
        return True


def _dumps_nodes(g):
    """Dump string of transition system states.

    @type g: L{FTS}
    """
    r = list()
    for u, d in g.nodes_iter(data=True):
        s = '\t Node: {u}, {values}\n'.format(
            u=u,
            values=', '.join(
                '{k} = {v}'.format(k=k, v=v)
                for k, v in d.iteritems()))
        r.append(s)
    return ''.join(r)
