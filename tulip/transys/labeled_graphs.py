# Copyright (c) 2013-2014 by California Institute of Technology
# and 2014 The Regents of the University of Michigan
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
# 3. Neither the name of the copyright holder(s) nor the names of its
#    contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
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
"""Rooted graphs and utils for labeled digraphs."""
from __future__ import absolute_import
import logging
logger = logging.getLogger(__name__)
import os
import networkx as nx


class SystemGraph(nx.MultiDiGraph):
    """Rooted multi-digraph with consistency methods.

    The attribute `initial_nodes` must contain nodes.
    """

    def __init__(self):
        super(SystemGraph, self).__init__()
        self.initial_nodes = set()

    def is_consistent(self):
        """Return `True` if attributes are conformant."""
        assert self.initial_nodes.issubset(self)

    def make_consistent(self):
        """Remove from attributes elements that are not nodes.

        Call this to update attributes like `initial_nodes`
        after removing nodes from `self`.
        """
        self.initial_nodes = self.initial_nodes.intersection(self)

    def dump(self, filename='out.pdf', rankdir='LR', prog='dot'):
        """Write to file using GraphViz.

        Requires
        ========
          - graphviz dot: http://www.graphviz.org/
          - pydot: https://pypi.python.org/pypi/pydot

        See Also
        ========
        `pydot.Dot.write`

        @param filename: path with exrension, for example: 'out.pdf'
        @type filename: str
        @param rankdir, prog: see `pydot.Dot.write`
        @param wrap: max width of node strings
        @type wrap: int

        @rtype: bool
        @return: True if saving completed successfully, False otherwise.
        """
        name, ext = os.path.splitext(filename)
        pd = self.to_pydot()
        pd.set_rankdir(rankdir)
        pd.set_splines('true')
        pd.write(filename, format=ext[1:], prog=prog)
        return True

    def to_pydot(self):
        return nx.to_pydot(self)

    def plot(self, rankdir='LR', prog=None, wrap=10, ax=None):
        """Plot image using GraphViz.

        No file I/O involved.
        Requires GraphViz dot and either Matplotlib or IPython.

        NetworkX does not yet support plotting multiple edges between 2 nodes.
        This method fixes that issue, so users don't need to look at files
        in a separate viewer during development.

        See Also
        ========
        L{dump}

        Depends
        =======
        dot and either of IPython or Matplotlib
        """
        # anything to plot ?
        if not self:
            print(
                60 * '!' +
                "\nThe system doesn't have any states to plot.\n"
                + 60 * '!')
            return
        if prog is None:
            prog = self.default_layout
        from tulip.transys.export import graph2dot
        return graph2dot.plot_pydot(self, prog, rankdir, wrap, ax=ax)


def find_nodes(self, states=None, with_attr_dict=None, **with_attr):
    """Filter by desired states and by desired state labels.

    Examples
    ========
    Assume that the system is:

    >>> import transys as trs
    >>> ts = trs.FTS()
    >>> ts.atomic_propositions.add('p')
    >>> ts.add_node('s0', ap={'p'})

      - To find all states with a specific label C{{'p'}}:

          >>> ts.add_node('s1', ap={'p'})
          >>> b = find_nodes(with_attr_dict={'ap':{'p'} } )
          >>> states = [state for (state, label_) in b]
          >>> print(set(states) )
          {'s0', 's1'}

      - To find all states in subset C{M} labeled with C{{'p'}}:

          >>> ts.add_node('s2', ap={'p'})
          >>> M = {'s0', 's2'}
          >>> b = find_nodes(M, {'ap': {'p'} } )
          >>> states = [state for (state, label_) in b]
          >>> print(set(states) )
          {'s0', 's2'}

    @param states: subset of states over which to search
    @type states: 'any' (default)
        | iterable of valid states
        | single valid state

    @param with_attr_dict: label with which to filter the states
    @type with_attr_dict: {sublabel_type : desired_sublabel_value, ...}
        | leave empty, to allow any state label (default)

    @param with_attr: label key-value pairs which take
        precedence over C{with_attr_dict}.

    @rtype: list of labeled states
    @return: [(C{state}, C{label}),...]
        where:
            - C{state} \\in C{states}
            - C{label}: dict
    """
    if with_attr_dict is None:
        with_attr_dict = with_attr
    else:
        try:
            with_attr_dict.update(with_attr)
        except AttributeError:
            raise Exception('with_attr_dict must be a dict')
    if states is not None:
        # singleton check
        if states in self:
            state = states
            msg = (
                'LabeledStates.find got single state: ' +
                str(state) + '\n'
                'instead of Iterable of states.\n')
            states = [state]
            msg += 'Replaced given states = ' + str(state)
            msg += ' with states = ' + str(states)
            logger.debug(msg)
    found_state_label_pairs = []
    for state, attr_dict in self.graph.nodes_iter(data=True):
        logger.debug('Checking state_id = ' + str(state) +
                     ', with attr_dict = ' + str(attr_dict))
        if states is not None:
            if state not in states:
                logger.debug('state_id = ' + str(state) + ', not desired.')
                continue
        msg = (
            'Checking state label:\n\t attr_dict = ' +
            str(attr_dict) +
            '\n vs:\n\t desired_label = ' + str(with_attr_dict))
        logger.debug(msg)
        if not with_attr_dict:
            logger.debug('Any label acceptable.')
            ok = True
        else:
            ok = label_is_desired(attr_dict, with_attr_dict)
        if ok:
            logger.debug('Label Matched:\n\t' + str(attr_dict) +
                         ' == ' + str(with_attr_dict))
            state_label_pair = (state, dict(attr_dict))
            found_state_label_pairs.append(state_label_pair)
        else:
            logger.debug('No match for label---> state discarded.')
    return found_state_label_pairs


def find_edges(self, from_states=None, to_states=None,
               with_attr_dict=None, typed_only=False, **with_attr):
    """Find all edges between given states with given labels.

    Instead of having two separate methods to:

      - find all labels of edges between given states (s1, s2)

      - find all transitions (s1, s2, L) with given label L,
            possibly from some given state s1,
            i.e., the edges leading to the successor states
            Post(s1, a) = Post(s1) restricted by action a

    this method provides both functionalities.

    Preimage under edge labeling function L of given label,
    intersected with given subset of edges::
        L^{-1}(desired_label) \\cap (from_states x to_states)

    See Also
    ========
    L{add}, L{add_adj}

    @param from_states: edges must start from this subset of states
    @type from_states:
        - iterable of existing states, or
        - None (no constraint, default)

    @param to_states: edges must end in this subset of states
    @type to_states:
        - iterable of existing states, or
        - None (no constraint, default)

    @param with_attr_dict: edges must be annotated with these labels
    @type with_attr_dict:
        - {label_type : desired_label_value, ...}, or
        - None (no constraint, default)

    @param with_attr: label type-value pairs,
        take precedence over C{desired_label}.

    @return: set of transitions = labeled edges::
            (from_state, to_state, label)
    such that::
            (from_state, to_state )
            in from_states x to_states

    @rtype: list of transitions::
            = list of labeled edges
            = [(from_state, to_state, label),...]
    where:
      - C{from_state} in C{from_states}
      - C{to_state} in C{to_states}
      - C{label}: dict
    """
    if with_attr_dict is None:
        with_attr_dict = with_attr
    try:
        with_attr_dict.update(with_attr)
    except:
        raise TypeError('with_attr_dict must be a dict')
    found_transitions = []
    u_v_edges = self.graph.edges_iter(nbunch=from_states, data=True)
    if to_states is not None:
        u_v_edges = [(u, v, d)
                     for u, v, d in u_v_edges
                     if v in to_states]
    for u, v, attr_dict in u_v_edges:
        ok = True
        if not with_attr_dict:
            logger.debug('Any label is allowed.')
        elif not attr_dict:
            logger.debug('No labels defined.')
        else:
            logger.debug('Checking guard.')
            ok = label_is_desired(attr_dict, with_attr_dict)
        if ok:
            logger.debug('Transition label matched desired label.')
            transition = (u, v, dict(attr_dict))
            found_transitions.append(transition)
    return found_transitions


def remove_deadends(g):
    """Recursively delete nodes with no outgoing transitions."""
    s = {1}
    while s:
        s = {n for n in g if not g.succ[n]}
        g.remove_nodes_from(s)
    g.make_consistent()
    return g


def paint(g, u, color):
    """Color the given state.

    The state is filled with given color,
    rendered with dot when plotting and saving.

    @param state: valid system state

    @param color: with which to paint C{state}
    @type color: str of valid dot color
    """
    g.add_node(u, style='filled', fillcolor=color)


def add_adj(g, adj, adj2states, **attr):
    """Add multiple labeled transitions from adjacency matrix.

    @param adj: new edges as adjacency matrix
    @type adj: `scipy.sparse.lil` (list of lists)
    @param adj2states: map from `adj` indices to `g` nodes.
    @type adj2states: `dict`
    """
    # square ?
    if adj.shape[0] != adj.shape[1]:
        raise Exception('Adjacency matrix must be square.')
    # convert to format friendly for edge iteration
    nx_adj = nx.from_scipy_sparse_matrix(adj, create_using=nx.DiGraph())
    # add each edge
    for i, j in nx_adj.edges_iter():
        u = adj2states[i]
        v = adj2states[j]
        g.add_edge(u, v, **attr)


def check_value(v, dom):
    if isinstance(dom, (set, list)):
        assert v in dom
    elif dom in {'bool', 'boolean'}:
        assert v in {0, 1, False, True}
    elif isinstance(dom, tuple):
        assert dom[0] <= int(v) <= dom[1]
    else:
        raise TypeError('dom not a `set` or `tuple` or "bool".')
