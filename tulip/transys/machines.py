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
"""Finite State Machines Module"""
from __future__ import absolute_import
import copy
from pprint import pformat
from random import choice
import networkx as nx
from tulip.transys.labeled_graphs import SystemGraph

# inline imports:
#
# from tulip.transys.export import machine2scxml

_hl = 40 * '-'

# port type
pure = {'present', 'absent'}


def is_valuation(ports, valuations):
    for name, port_type in ports.iteritems():
        curvaluation = valuations[name]
        # functional set membership description ?
        if callable(port_type):
            ok = port_type(curvaluation)
        else:
            ok = curvaluation in port_type
        if not ok:
            raise TypeError('Not a valuation.')


def create_machine_ports(spc_vars):
    """Create proper port domains of valuations, given port types.

    @param spc_vars: port names and types inside tulip.
        For arbitrary finite types the type can be a list of strings,
        instead of a range of integers.
        These are as originally defined by the user or synth.
    """
    ports = dict()
    for env_var, var_type in spc_vars.iteritems():
        if var_type == 'boolean':
            domain = {0, 1}
        elif isinstance(var_type, tuple):
            # integer domain
            start, end = var_type
            domain = set(range(start, end + 1))
        elif isinstance(var_type, list):
            # arbitrary finite domain defined by list var_type
            domain = set(var_type)
        ports[env_var] = domain
    return ports


class Transducer(SystemGraph):
    """Sequential Transducer, i.e., a letter-to-letter function.

    Inputs
    ======
    P = {p1, p2,...} is the set of input ports.
    An input port p takes values in a set Vp.
    Set Vp is called the "type" of input port p.
    A "valuation" is an assignment of values to the input ports in P.

    We call "inputs" the set of pairs::

      {(p_i, Vp_i),...}

    of input ports p_i and their corresponding types Vp_i.

    A guard is a predicate (bool-valued) used as sub-label for a transition.
    A guard is defined by a set and evaluated using set membership.
    So given an input port value p=x, then if::

      x \in guard_set

    then the guard is True, otherwise it is False.

    The "inputs" are defined by an OrderedDict::

      {'p1':explicit, 'p2':check, 'p3':None, ...}

    where:
      - C{explicit}:
        is an iterable representation of Vp,
        possible only for discrete Vp.
        If 'p1' is explicitly typed, then guards are evaluated directly::

          input_port_value == guard_value ?

      - C{check}:
        is a class with methods:

          - C{__contains__(x) }:
            check if guard value given to input port 'p1' is
            in the set of possible values Vp.

          - C{__call__(guard_set, input_port_value) }:
            check if C{input_port_value} \\in C{guard_set}
            This allows symbolic type definitions.

            For example, C{input_port_value} might be assigned
            int values, but the C{guard_set} be defined by
            a symbolic expression as the str: 'x<=5'.

            Then the user is responsible for providing
            the appropriate method to the Mealy Machine,
            using the custom C{check} class described here.

            Note that we could provide a rudimentary library
            for the basic types of checks, e.g., for
            the above simple symbolic case, where using
            function eval() is sufficient.

      - C{None}:
        signifies that no type is currently defined for
        this input port, so input type checking and guard
        evaluation are disabled.

        This can be used to skip type definitions when
        they are not needed by the user.

        However, since Machines are in general the output
        of synthesis, it follows that they are constructed
        by code, so the benefits of typedefs will be
        considerable compared to the required coding effort.

    Guards annotate edges::

      Guards: States x States ---> Input_Predicates

    Outputs
    =======
    Similarly defined to inputs, but:

      - for Mealy Machines they annotate edges
      - for Moore Machines they annotate nodes

    State Variables
    ===============
    Similarly defined to inputs, they annotate states,
    for both Mealy and Moore machines::

      States ---> State_Variables

    Update Function
    ===============
    The transition relation:

      - for Mealy Machines::

        States x Input_Valuations ---> Output_Valuations x States

        Note that in the range Output_Valuations are ordered before States
        to emphasize that an output_valuation is produced
        during the transition, NOT at the next state.

        The data structure representation of the update function is
        by storage of the Guards function and definition of Guard
        evaluation for each input port via the OrderedDict discussed above.

      - for Moore Machines::

        States x Input_Valuations ---> States
        States ---> Output_valuations

    Note
    ====
    A transducer may operate on either finite or infinite words, i.e.,
    it is not equipped with interpretation semantics on the words,
    so it does not "care" about word length.
    It continues as long as its input is fed with letters.

    For Machines, each state label consists of (possibly multiple) sublabels,
    each of which is either a variable, or, only for Moore machines,
    may be an output.

    See Also
    ========
    FSM, MealyMachine, MooreMachine
    """

    def __init__(self):
        # values will point to values of _*_label_def below
        self.state_vars = dict()
        self.inputs = dict()
        self.outputs = dict()
        # self.set_actions = {}
        super(Transducer, self).__init__()


class MooreMachine(Transducer):
    """Moore machine.

    A Moore machine implements the discrete dynamics::
        x[k+1] = f(x[k], u[k] )
        y[k] = g(x[k] )
    where:
      - k: discrete time = sequence index
      - x: state = valuation of state variables
      - X: set of states = S
      - u: inputs = valuation of input ports
      - y: output actions = valuation of output ports
      - f: X-> 2^X, transition function
      - g: X-> Out, output function
    Observe that the output depends only on the state.

    Note
    ====
    valuation: assignment of values to each port

    Reference
    =========
    U{[M56]
    <http://tulip-control.sourceforge.net/doc/bibliography.html#m56>}
    """

    def __str__(self):
        """Get informal string representation."""
        s = (
            _hl + '\nMoore Machine: ' + self.name + '\n' + _hl + '\n' +
            'State Variables:\n\t(name : type)\n' +
            _print_ports(self.state_vars) +
            'Input Ports:\n\t(name : type)\n' +
            _print_ports(self.inputs) +
            'Output Ports:\n\t(name : type)\n' +
            _print_ports(self.outputs) +
            'States & State Var Values: (state : outputs : vars)\n')
        for state, label_dict in self.nodes_iter(data=True):
            s += '\t' + str(state) + ' :\n'
            # split into vars and outputs
            var_values = {k: v for k, v in label_dict.iteritems()
                          if k in self.state_vars}
            output_values = {k: v for k, v in label_dict.iteritems()
                             if k in self.outputs}
            s += (_print_label(var_values) + ' : ' +
                  _print_label(output_values))
        s += (
            'Initial States:\n' +
            pformat(self.initial_nodes, indent=3) + 2 * '\n')
        s += 'Edges & Labels: (from --> to : label)\n'
        for u, v, d in self.edges_iter(data=True):
            s += (
                '\t' + str(u) + ' ---> ' +
                str(v) + ' :\n' +
                _print_label(d))
        s += _hl + '\n'
        return s

    def to_pydot(self):
        g = nx.MultiDiGraph()
        for u, d in self.nodes_iter(data=True):
            var = _join(d, self.state_vars)
            o = _join(d, self.outputs)
            r = list()
            if var:
                r.append('var:\n' + var)
            if o:
                r.append('out:\n' + o)
            label = '\n'.join(r)
            g.add_node(u, label=label, shape='ellipse')
        for u, v, d in self.edges_iter(data=True):
            i = _join(d, self.inputs)
            if i:
                label = 'in:\n' + i
            else:
                label = ''
            g.add_edge(u, v, label=label)
        return nx.to_pydot(g)


class MealyMachine(Transducer):
    """Mealy machine.

    Examples
    ========
    Traffic Light: Fig. 3.14, p.72 U{[LS11]
    <http://tulip-control.sourceforge.net/doc/bibliography.html#ls11>}

    >>> m = MealyMachine()
    >>> pure_signal = {'present', 'absent'}
    >>> m.inputs.update([('tick', pure_signal) ])
    >>> m.outputs.update([('go', pure_signal), ('stop', pure_signal) ])
    >>> m.add_nodes_from(['red', 'green', 'yellow'])
    >>> m.initial_nodes.add('red')

    For brevity:

    >>> p = 'present'
    >>> a = 'absent'

    >>> label = {'tick':p, 'go':p, 'stop':a}
    >>> m.add_edge('red', 'green', **label)
    >>> label = {'tick':p, 'go':a, 'stop':p}
    >>> m.add_edge('green', 'yellow', **label)
    >>> label = {'tick':p, 'go':a, 'stop':p}
    >>> m.add_edge('yellow', 'red', **label)

    Theory
    ======
    A Mealy machine implements the discrete dynamics::
        x[k+1] = f(x[k], u[k] )
        y[k] = g(x[k], u[k] )
    where:
      - k: discrete time = sequence index
      - x: state = valuation of state variables
      - X: set of states = S
      - u: inputs = valuation of input ports
      - y: output actions = valuation of output ports
      - f: X-> 2^X, transition function
      - g: X-> Out, output function
    Observe that the output is defined when a reaction occurs to an input.

    Note
    ====
    valuation: assignment of values to each port

    Reference
    =========
    U{[M55]
    <http://tulip-control.sourceforge.net/doc/bibliography.html#m55>}
    """

    def __str__(self):
        """Get informal string representation."""
        s = (
            _hl + '\nMealy Machine: ' + self.name + '\n' + _hl + '\n' +
            'State Variables:\n\t(name : type)\n' +
            _print_ports(self.state_vars))
        s += 'States & State Var Values:\n'
        for state, label_dict in self.nodes_iter(data=True):
            s += ('\t' + str(state) + ' :\n' +
                  _print_label(label_dict))
        s += (
            'Initial States:\n' +
            pformat(self.initial_nodes, indent=3) + 2 * '\n' +
            'Input Ports:\n\t(name : type)\n' +
            _print_ports(self.inputs) +
            'Output Ports:\n\t(name : type)\n' +
            _print_ports(self.outputs) +
            'Edges & Labels: (from --> to : label)\n')
        for u, v, d in self.edges_iter(data=True):
            s += (
                '\t' + str(u) + ' ---> ' +
                str(v) + ' :\n' +
                _print_label(d))
        s += _hl + '\n'
        return s

    def to_pydot(self, ports=None):
        g = nx.MultiDiGraph()
        for u, d in self.nodes_iter(data=True):
            label = '{u}\n{l}'.format(u=u, l=_join(d, self.state_vars))
            g.add_node(u, label=label, shape='ellipse')
        for u, v, d in self.edges_iter(data=True):
            if ports is None:
                din = self.inputs
                dout = self.outputs
            else:
                din = project_dict(self.inputs, ports)
                dout = project_dict(self.outputs, ports)
            i = _join(d, din)
            o = _join(d, dout)
            r = list()
            if i:
                r.append('in:\n' + i)
            if o:
                r.append('out:\n' + o)
            label = '\n'.join(r)
            label = label if label else '""'
            g.add_edge(u, v, label=label)
        return nx.to_pydot(g)

    def reaction(self, from_state, inputs):
        """Return next state and output, when reacting to given inputs.

        The machine must be deterministic.
        (for each state and input at most a single transition enabled,
        this notion does not coincide with output-determinism)

        Not exactly a wrapper of L{find_edges},
        because it matches only that part of an edge label
        that corresponds to the inputs.

        @param from_state: source node of transition

        @param inputs: C{dict} assigning a valid value to each input port.
        @type inputs: {'port_name': port_value, ...}

        @return: output values and next state.
        @rtype: (outputs, next_state)
          where C{outputs}: C{{'port_name': port_value, ...}}
        """
        # match only inputs (explicit valuations, not symbolic)
        enabled_trans = [
            (i, j, d)
            for i, j, d in self.edges_iter([from_state], data=True)
            if project_dict(d, self.inputs) == inputs]
        # must be deterministic
        try:
            ((_, next_state, attr_dict), ) = enabled_trans
        except ValueError:
            raise Exception(
                'must be input-deterministic, '
                'found enabled transitions: '
                '{t}'.format(t=enabled_trans))
        outputs = project_dict(attr_dict, self.outputs)
        return (next_state, outputs)

    def run(self, from_state=None, input_sequences=None):
        """Guided or interactive run.

        @param input_sequences: if C{None}, then call L{interactive_run},
            otherwise call L{guided_run}.

        @return: output of L{guided_run}, otherwise C{None}.
        """
        if input_sequences is None:
            interactive_run(self, from_state=from_state)
        else:
            return guided_run(self, from_state=from_state,
                              input_sequences=input_sequences)


# note on non-determinism and simulation:
# =====
#
# Transducers here have deterministic semantics,
# even if the exist non-deterministic choices of
# After all, this turns out to be the essentialy distinction
# between a non-deterministic game graph (generator) and
# a machine (transducer).
# (Modulo some graph and labeling transformations.)
#
# A Game Graph (or a transition system if only one player exists,
# so no partition on the states accompanies the labeled directed graph)
# is a generator, which takes "all possible transitions".
# Same for acceptors, in which case the resulting language is filtered
# by the accepting condition, to yield the language
# represented by the acceptor.
#
# In other words: syntactically they are more or less the same,
# but semantically they are different.
# Therefore the methods or functions that manipulate them differ,
# because they involve their semantics.
#
# A possible alternative in the future would be to
#   - raise an exception in case of non-determinism
#   - unless the user explicitly allows non-determinism
# In the latter case, machine non-determinism is resolved
# arbitrarily during simulation.
# I think this would correspond to what happens to random simulation
# in SPIN also.

# note on impossibility to simulate non-deterministic machine
# =====
#
# It is currently impossible to simulate a non-deterministic machine.
# A non-deterministic machine represents a set of transducers that
# are all valid solutions to the synthesis problem.
#
# But because the semantics of real-world execution require that only
# a single trace be produced (so resolve the non-determinism in an
# arbitrary way), this use is discouraged.
#
# In the future, such an option could be added.
# However, it may be a bad design approach.
# Instead, since a non-deterministic machine is a set of satisfying models
# (= solutions), a better approach would be to first fix one deterministic
# alternative as the desired solution, then simulate it.
#
# Finally there also exists the issue of initial condition interpretation.
# Formally, if the initial condition of sys is interpreted as an env var,
# then after we learn the initial system state of the particular instance
# that will be simulated, a deterministic transducer should be obtained
# with only that state as initial, then simulated.
#
# A possible alternative in the future would be to
# allow for initial system state
# non-determinism, but require that the initial system state be passed as
# an argument in this case.


# note on dict of lists vs list of dicts
# =====
#
# Typically there are few ports but many time indices
# so checking for equal len of all input sequences is cheaper
# than checking no port is missing for each dict in a list of dicts over time.
#
# Also less memory is needed (one dict of a few lists of 1000 items each,
#   not 1000 dicts with the same keys of a few items each)
#
# Absent inputs at certain times can still be represented,
# by explicictly assigning them the value None.
# This avoids mistakes (explicit better than inplicit, cf Zen of Python)
# and ensures uniformity in the data structure (no absent terms)
#
# Also, the dict of lists representation is friendly for adding/removing
# ports. This may prove convenient if working with a network of machines.
# For example splitting many output sequences
# to become inputs for other machines.

# note on terminology
# =====
#
# The term "simulation" is ambiguous and not used any more:
#
# 1. it is too general
# 2. it is dangerous, because it has also another meaning in this context
# 3. "run" is also wrongly used here (i.e., more than the sequence )


def guided_run(mealy, from_state=None, input_sequences=None):
    """Run deterministic machine reacting to given inputs.

    @param from_state: start simulation

    @param mealy: input-deterministic Mealy machine
    @type mealy: L{MealyMachine}

    @param from_state: start simulation at this state.
        If C{None}, then use the unique initial state C{Sinit}.

    @param input_sequences: one sequence of values for each input port
    @type input_sequences: C{dict} of C{lists}

    @return: sequence of states and sequence of output valuations
    @rtype: (states, output_sequences)
      where:
        - C{states} is a C{list} of states excluding C{from_state}
        - C{output_sequences} is a C{dict} of C{lists}
    """
    seqs = input_sequences  # abbrv
    missing_ports = set(mealy.inputs).difference(seqs)
    if missing_ports:
        raise ValueError('missing input port(s): ' + missing_ports)
    # dict of lists ?
    non_lists = {k: v for k, v in seqs.iteritems() if not isinstance(v, list)}
    if non_lists:
        raise TypeError('Values must be lists, for: ' + str(non_lists))
    # uniform list len ?
    if len(set(len(x) for x in seqs.itervalues())) > 1:
        raise ValueError('All input sequences must be of equal length.')
    # note: initial sys state non-determinism not checked
    # initial sys edge non-determinism checked instead (more restrictive)
    if from_state is None:
        state = next(iter(mealy.initial_nodes))
    else:
        state = from_state
    n = len(next(seqs.itervalues()))
    states_seq = []
    output_seqs = {k: list() for k in mealy.outputs}
    for i in range(n):
        inputs = {k: v[i] for k, v in seqs.iteritems()}
        state, outputs = mealy.reaction(state, inputs)
        states_seq.append(state)
        for k in output_seqs:
            output_seqs[k].append(outputs[k])
    return (states_seq, output_seqs)


def random_run(mealy, from_state=None, N=10):
    """Return run from given state for N random inputs.

    Inputs selected randomly in a way that does not block the machine
    So they are not arbitrarily random.
    If the machine is a valid synthesis solution,
    then all safe environment inputs can be generated this way.

    Randomly generated inputs may violate liveness assumption on environment.

    @param mealy: input-deterministic Mealy machine
    @type mealy: C{MealyMachine}

    @param N: number of reactions (inputs)
    @type N: int

    @return: same as L{guided_run}
    """
    if from_state is None:
        state = next(iter(mealy.initial_nodes))
    else:
        state = from_state
    states_seq = []
    output_seqs = {k: list() for k in mealy.outputs}
    for i in xrange(N):
        trans = mealy.out_edges(state)
        # choose next transition
        selected_trans = choice(trans)
        _, new_state, attr_dict = selected_trans
        # extend execution trace
        states_seq.append(new_state)
        # extend output traces
        outputs = project_dict(attr_dict, mealy.outputs)
        for k in output_seqs:
            output_seqs[k].append(outputs[k])
        # updates
        old_state = state
        state = new_state
        # printing
        inputs = project_dict(attr_dict, mealy.inputs)
        print(
            'move from\n\t state: ' + str(old_state) +
            '\n\t with input:' + str(inputs) +
            '\n\t to state: ' + str(new_state) +
            '\n\t reacting by producing output: ' + str(outputs))
    return (states_seq, output_seqs)


def interactive_run(mealy, from_state=None):
    """Run input-deterministic Mealy machine using user input.

    @param mealy: input-deterministic Mealy machine
    @type mealy: L{MealyMachine}
    """
    if from_state is None:
        state = next(iter(mealy.initial_nodes))
    else:
        state = from_state
    while True:
        print('\n Current state: ' + str(state))
        if _interactive_run_step(mealy) is None:
            break


def _interactive_run_step(mealy, state):
    if state is None:
        raise Exception('Current state is None')
    # note: the spaghettiness of previous version was caused
    #   by interactive simulation allowing both output-non-determinism
    #   and implementing spawning (which makes sense only for generators,
    #   *not* for transducers)
    trans = mealy.out_edges(state)
    if not trans:
        print('Stop: no outgoing transitions.')
        return None
    while True:
        try:
            selected_trans = _select_transition(mealy, trans)
        except:
            print('Selection not recognized. Please try again.')
    if selected_trans is None:
        return None
    (from_, to_state, attr_dict) = selected_trans
    inputs = project_dict(attr_dict, mealy.inputs)
    outputs = project_dict(attr_dict, mealy.outputs)
    print(
        'Moving from state: ' + str(state) +
        ', to state: ' + str(to_state) + '\n' +
        'given inputs: ' + str(inputs) + '\n' +
        'reacting with outputs: ' + str(outputs))
    return True


def _select_transition(mealy, trans):
    msg = 'Found more than 1 outgoing transitions:' + 2 * '\n'
    for i, t in enumerate(trans):
        (from_state, to_state, attr_dict) = t
        inputs = project_dict(attr_dict, mealy.inputs)
        outputs = project_dict(attr_dict, mealy.outputs)
        msg += (
            '\t' + str(i) + ' : ' +
            str(from_state) + ' ---> ' + str(to_state) + '\n' +
            '\t inputs:' + str(inputs) +
            '\t outputs:' + str(outputs) +
            '\n\n')
    msg += (
        '\n' +
        'Select from the available transitions above\n' +
        'by giving its integer,\n' +
        'Press "Enter" to stop the simulation:\n' +
        '\t int = ')
    id_selected = raw_input(msg)
    if not id_selected:
        return None
    return trans[int(id_selected)]


def moore2mealy(moore):
    """Convert Moore machine to equivalent Mealy machine.

    Reference
    =========
    U{[LS11]
    <http://tulip-control.sourceforge.net/doc/bibliography.html#ls11>}

    @type moore: L{MooreMachine}

    @rtype: L{MealyMachine}
    """
    if not isinstance(moore, MooreMachine):
        raise TypeError('moore must be a MooreMachine')
    mealy = MealyMachine()
    mealy.inputs.update(mealy.inputs)
    mealy.outputs.update(moore.outputs)
    mealy.add_nodes_from(moore)
    mealy.initial_nodes.update(moore.initial_nodes)
    # cp transitions
    for si in moore:
        output_values = {
            k: v for k, v in moore[si].iteritems()
            if k in moore.outputs}
        output_values = copy.deepcopy(output_values)
        for si_, sj, d in moore.out_edges_iter(si):
            # note that we don't filter only input ports,
            # so other edge annotation is preserved
            q = dict(d)
            q.update(output_values)
            mealy.add_edge(si, sj, q)
    return mealy


def mealy2moore(mealy):
    """Convert Mealy machine to almost equivalent Moore machine.

    A Mealy machine cannot be transformed to an equivalent Moore machine.
    It can be converted to a Moore machine with an arbitrary initial output,
    which outputs the Mealy output at its next reaction.

    Reference
    =========
    U{[LS11]
    <http://tulip-control.sourceforge.net/doc/bibliography.html#ls11>}

    @type mealy: L{MealyMachine}

    @rtype: L{MooreMachine}
    """
    # TODO: check for when Mealy is exactly convertible to Moore
    if not isinstance(mealy, MealyMachine):
        raise TypeError('moore must be a MealyMachine')
    moore = MooreMachine()
    moore.inputs.update(mealy.inputs)
    moore.outputs.update(mealy.outputs)
    # initial state with arbitrary label
    out = {k: list(v)[0] for k, v in mealy.outputs.iteritems()}
    s0 = next(iter(mealy.initial_nodes))
    # create maps between Moore and Mealy states
    moore2mealy_states = dict()  # {qj : si} (function)
    mealy2moore_states = dict()  # {si : {qj, qk, ...} } (relation)
    new_s0 = _create_state_str(
        s0, out, moore, moore2mealy_states,
        mealy2moore_states)
    moore.add_node(new_s0, out)
    moore.initial_nodes.add(new_s0)
    # cp transitions and create appropriate states
    Q = set()
    S = set()
    Q.add(new_s0)
    S.add(new_s0)
    while Q:
        u = Q.pop()
        si = moore2mealy_states[u]
        for _, oldv, d in mealy.edges_iter(si, data=True):
            in_values, out_values = _split_io(d, mealy)
            v = _create_state_str(
                oldv, out_values, moore, moore2mealy_states,
                mealy2moore_states)
            moore.add_edge(u, v, in_values)
            if v not in S:
                Q.add(v)
                S.add(v)
    return moore


def _print_ports(port_dict):
    s = ''
    for port_name, port_type in port_dict.iteritems():
        s += '\t' + str(port_name) + ' : '
        s += pformat(port_type) + '\n'
    s += '\n'
    return s


def _print_label(label_dict):
    s = ''
    for name, value in label_dict.iteritems():
        s += '\t\t' + str(name) + ' : ' + str(value) + '\n'
    s += '\n'
    return s


def _create_state_str(mealy_state, output, moore,
                      moore2mealy_states,
                      mealy2moore_states):
    """Used to create Moore states when converting Mealy -> Moore."""
    for s in mealy2moore_states.setdefault(mealy_state, set()):
        # check output values
        if moore.node[s] == output:
            return s
    # create new
    n = len(moore)
    s = 's' + str(n)
    moore.add_node(s, output)
    moore2mealy_states[s] = mealy_state
    mealy2moore_states[mealy_state].add(s)
    return s


def _split_io(attr_dict, machine):
    """Split into inputs and outputs."""
    input_values = {k: v for k, v in attr_dict.iteritems()
                    if k in machine.inputs}
    output_values = {k: v for k, v in attr_dict.iteritems()
                     if k in machine.outputs}
    return input_values, output_values


project_dict = lambda x, y: {k: x[k] for k in x if k in y}
trim_dict = lambda x, y: {k: x[k] for k in x if k not in y}


def strip_ports(mealy, names):
    """Remove ports in C{names}.

    For example, to remove the atomic propositions
    labeling the transition system C{ts} used
    (so they are dependent variables), call it as:

      >>> strip_ports(mealy, ts.atomic_propositions)

    @type mealy: L{MealyMachine}

    @type names: iterable container of C{str}
    """
    new = MealyMachine()

    new.inputs.update(trim_dict(mealy.inputs, names))
    new.outputs.update(trim_dict(mealy.outputs, names))

    new.add_nodes_from(mealy)
    new.initial_nodes.update(mealy.initial_nodes)

    for u, v, d in mealy.edges_iter(data=True):
        d = trim_dict(d, names)
        new.add_edge(u, v, **d)
    return new


def _join(d, keys, sep=': ', itemsep='\n'):
    return itemsep.join(
        '{k}{sep}{v}'.format(k=k, sep=sep, v=v)
        for k, v in d.iteritems()
        if k in keys)
