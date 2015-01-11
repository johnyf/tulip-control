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
"""Base classes for labeled directed graphs"""
from __future__ import absolute_import
import logging
logger = logging.getLogger(__name__)
import os
import copy
from pprint import pformat
from collections import Iterable
import warnings
import networkx as nx
from tulip.transys.mathset import SubSet, TypedDict

# inline imports:
#
# from tulip.transys.export import graph2dot
# from tulip.transys.export import save_d3
# from tulip.transys.export import graph2dot




    def paint(self, state, color):
        """Color the given state.

        The state is filled with given color,
        rendered with dot when plotting and saving.

        @param state: valid system state

        @param color: with which to paint C{state}
        @type color: str of valid dot color
        """
        self.graph.node[state]['style'] = 'filled'
        self.graph.node[state]['fillcolor'] = color

    def find(self, states=None, with_attr_dict=None, **with_attr):
        """Filter by desired states and by desired state labels.

        Examples
        ========
        Assume that the system is:

        >>> import transys as trs
        >>> ts = trs.FTS()
        >>> ts.atomic_propositions.add('p')
        >>> ts.states.add('s0', ap={'p'})

          - To find the label of a single state C{'s0'}:

              >>> a = ts.states.find(['s0'] )
              >>> (s0_, label) = a[0]
              >>> print(label)
              {'ap': set(['p'])}

          - To find all states with a specific label C{{'p'}}:

              >>> ts.states.add('s1', ap={'p'})
              >>> b = ts.states.find(with_attr_dict={'ap':{'p'} } )
              >>> states = [state for (state, label_) in b]
              >>> print(set(states) )
              {'s0', 's1'}

          - To find all states in subset C{M} labeled with C{{'p'}}:

              >>> ts.states.add('s2', ap={'p'})
              >>> M = {'s0', 's2'}
              >>> b = ts.states.find(M, {'ap': {'p'} } )
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

    def add_adj(
            self, adj, adj2states, attr_dict=None,
            check=True, **attr
    ):
        """Add multiple labeled transitions from adjacency matrix.

        The label can be empty.
        For more details see L{add}.

        @param adj: new transitions represented by adjacency matrix.
        @type adj: scipy.sparse.lil (list of lists)

        @param adj2states: map from adjacency matrix indices to states.
            If value not a state, raise Exception.
            Use L{States.add}, L{States.add_from} to add states first.

            For example the 1st state in adj2states corresponds to
            the first node in C{adj}.

            States must have been added using:

               - sys.states.add, or
               - sys.states.add_from

            If C{adj2states} includes a state not in sys.states,
            no transition is added and an exception raised.
        @type adj2states: either of:
            - C{dict} from adjacency matrix indices to
              existing, or
            - C{list} of existing states
        """
        # square ?
        if adj.shape[0] != adj.shape[1]:
            raise Exception('Adjacency matrix must be square.')
        # check states exist, before adding any transitions
        for state in adj2states:
            if state not in self.graph:
                raise Exception(
                    'State: ' + str(state) + ' not found.'
                    ' Consider adding it with sys.states.add')
        # convert to format friendly for edge iteration
        nx_adj = nx.from_scipy_sparse_matrix(
            adj, create_using=nx.DiGraph())
        # add each edge using existing checks
        for i, j in nx_adj.edges_iter():
            si = adj2states[i]
            sj = adj2states[j]
            self.add(si, sj, attr_dict, check, **attr)

    def find(self, from_states=None, to_states=None,
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
    def remove_deadends(self):
        """Recursively delete nodes with no outgoing transitions."""
        s = {1}
        while s:
            s = {n for n in self if not self.succ[n]}
            self.states.remove_from(s)

    def dot_str(self, wrap=10, **kwargs):
        """Return dot string.

        Requires pydot.
        """
        from tulip.transys.export import graph2dot
        return graph2dot.graph2dot_str(self, wrap, **kwargs)

    def save(self, filename=None, fileformat=None,
             rankdir='LR', prog=None,
             wrap=10, tikz=False):
        """Save image to file.

        Recommended file formats:

            - tikz (via dot2tex)
            - pdf
            - svg
            - dot
            - png

        Any other format supported by C{pydot.write} is available.

        Experimental:

            - html (uses d3.js)
            - 'scxml'

        Requires
        ========
          - graphviz dot: http://www.graphviz.org/
          - pydot: https://pypi.python.org/pypi/pydot

        and for tikz:

          - dot2tex: https://pypi.python.org/pypi/dot2tex
          - dot2texi: http://www.ctan.org/pkg/dot2texi
            (to automate inclusion)

        See Also
        ========
        L{plot}, C{pydot.Dot.write}

        @param filename: file path to save image to
            Default is C{self.name}, unless C{name} is empty,
            then use 'out.pdf'.

            If extension is missing '.pdf' is used.
        @type filename: str

        @param fileformat: replace the extension of C{filename}
            with this. For example::

                filename = 'fig.pdf'
                fileformat = 'svg'

            result in saving 'fig.svg'

        @param rankdir: direction for dot layout
        @type rankdir: str = 'TB' | 'LR'
            (i.e., Top->Bottom | Left->Right)

        @param prog: executable to call
        @type prog: dot | circo | ... see pydot.Dot.write

        @param wrap: max width of node strings
        @type wrap: int

        @param tikz: use tikz automata library in dot

        @rtype: bool
        @return: True if saving completed successfully, False otherwise.
        """
        if filename is None:
            if not self.name:
                filename = 'out'
            else:
                filename = self.name
        fname, fextension = os.path.splitext(filename)
        # default extension
        if not fextension or fextension is '.':
            fextension = '.pdf'
        if fileformat:
            fextension = '.' + fileformat
        filename = fname + fextension
        # drop '.'
        fileformat = fextension[1:]
        # check for html
        if fileformat is 'html':
            from tulip.transys.export import save_d3
            return save_d3.labeled_digraph2d3(self, filename)
        # subclass has extra export formats ?
        if hasattr(self, '_save'):
            if self._save(filename, fileformat):
                return True
        if prog is None:
            prog = self.default_layout
        from tulip.transys.export import graph2dot
        graph2dot.save_dot(self, filename, fileformat, rankdir,
                           prog, wrap, tikz=tikz)
        return True

    def plot(self, rankdir='LR', prog=None, wrap=10, ax=None):
        """Plot image using dot.

        No file I/O involved.
        Requires GraphViz dot and either Matplotlib or IPython.

        NetworkX does not yet support plotting multiple edges between 2 nodes.
        This method fixes that issue, so users don't need to look at files
        in a separate viewer during development.

        See Also
        ========
        L{save}

        Depends
        =======
        dot and either of IPython or Matplotlib
        """
        # anything to plot ?
        if not self.states:
            print(
                60 * '!' +
                "\nThe system doesn't have any states to plot.\n"
                + 60 * '!')
            return
        if prog is None:
            prog = self.default_layout
        from tulip.transys.export import graph2dot
        return graph2dot.plot_pydot(self, prog, rankdir, wrap, ax=ax)


def str2singleton(ap_label):
    """If string, convert to set(string).

    Convention: singleton str {'*'}
    can be passed as str '*' instead.
    """
    if isinstance(ap_label, str):
        logger.debug('Saw str state label:\n\t' + str(ap_label))
        ap_label = {ap_label}
        logger.debug('Replaced with singleton:\n\t' + str(ap_label) + '\n')
    return ap_label


def prepend_with(states, prepend_str):
    """Prepend items with given string.

    Example
    =======
    >>> states = [0, 1]
    >>> prepend_str = 's'
    >>> states = prepend_with(states, prepend_str)
    >>> assert(states == ['s0', 's1'] )

    See Also
    ========
    L{tuple2ba}, L{tuple2fts}

    @param states: items prepended with string C{prepend_str}
    @type states: iterable

    @param prepend_str: text prepended to C{states}.  If None, then
        C{states} is returned without modification
    @type prepend_str: str or None
    """
    if not isinstance(states, Iterable):
        raise TypeError('states must be Iterable. Got:\n\t' +
                        str(states) + '\ninstead.')
    if not isinstance(prepend_str, str) and prepend_str is not None:
        raise TypeError('prepend_str must be of type str. Got:\n\t' +
                        str(prepend_str) + '\ninstead.')
    if prepend_str is None:
        return states
    return [prepend_str + str(s) for s in states]
