# Copyright (c) 2012-2014 by California Institute of Technology
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
"""Interface to library of synthesis tools, e.g., JTLV, gr1c"""
from __future__ import absolute_import
import logging
logger = logging.getLogger(__name__)
import copy
import warnings
from tulip.transys import MealyMachine
from tulip.transys import machines
from tulip.transys.labeled_graphs import remove_deadends
from tulip.spec import GRSpec
from tulip.interfaces import jtlv, gr1c


_hl = '\n' + 60 * '-'


def _prime_dict(d):
    """Return `dict` with primed keys and `dict` mapping to them."""
    p = dict((_prime(k), d[k]) for k in d)
    d2p = {k: _prime(k) for k in d}
    return p, d2p


def _prime(s):
    return s + "'"


def _pstr(x):
    return '({x})'.format(x=x)


def _disj(s):
    return ' || '.join('({x})'.format(x=x) for x in s if x != '')


def _conj(s):
    return ' && '.join('({x})'.format(x=x) for x in s if x != '')


def _conj_intersection(s, r, paren=True):
    if paren:
        return ' && '.join('({x})'.format(x=x) for x in s if x in r)
    else:
        return ' && '.join('{x}'.format(x=x) for x in s if x in r)


def _conj_neg(s, paren=True):
    if paren:
        return ' && '.join('!({x})'.format(x=x) for x in s)
    else:
        return ' && '.join('!{x}'.format(x=x) for x in s)


def _conj_neg_diff(s, r, paren=True):
    if paren:
        return ' && '.join('({x})'.format(x=x) for x in s if x not in r)
    else:
        return ' && '.join('{x}'.format(x=x) for x in s if x not in r)


def mutex(iterable):
    """Mutual exclusion for all time."""
    iterable = filter(lambda x: x != '', iterable)
    if not iterable:
        return list()
    if len(iterable) <= 1:
        return []
    return [_conj([
        '!(' + str(x) + ') || (' + _conj_neg_diff(iterable, [x]) + ')'
        for x in iterable])]


def exactly_one(iterable):
    """N-ary xor.

    Contrast with pure mutual exclusion.
    """
    if len(iterable) <= 1:
        return [_pstr(x) for x in iterable]
    return ['(' + _disj([
        '(' + str(x) + ') && ' + _conj_neg_diff(iterable, [x])
        for x in iterable]) + ')']


# duplicate states are impossible, because each networkx node is unique
# non-contiguous integers for states fine: you are lossing efficiency
# -- synth doesn't care about that


def iter2var(var, values):
    """Assignments for finite domain.

    An integer or string variable can be used.

    If all values are integers, then an integer is used.
    If all values are strings, then a string variable is used.
    Otherwise an exception is raised, unless Booleans have been requested.

    @param values: domain of `int` or `str` variable.
    @type values: iterable container
    @param var: name to use for integer or string valued variabe.
    @type var: `str`

    @return: `tuple` of:
      - mapping from values to GR(1) actions.
        If Booleans are used, then GR(1) are the same.
        Otherwise, they map to e.g. 'act = "wait"' or 'act = 3'
    @rtype: `dict`
    """
    if not values:
        logger.debug('empty container, so empty dict for solver expr')
        return dict(), None
    logger.debug('mapping domain: {values}\n\t'.format(values=values))
    all_int = all(isinstance(x, int) for x in values)
    all_str = all(isinstance(x, str) for x in values)
    if all_int:
        domain = (min(values), max(values))
    elif all_str:
        domain = list(values)
    else:
        raise TypeError('Integer and string states must not be mixed.')
    return domain


def sys_to_spec(g, nodevar, ignore_initial, receptive=False):
    """Convert transition system to GR(1) fragment of LTL.

    The attribute `g.owner` defines who selects the next node.
    The attribute `g.env_vars` determines who controls each variable.

    Caution
    =======
    Initial values of variables that label only edges
    have to be specified separately.

    @param g: `TransitionSystem`
    @param nodevar: variable that stores current node
    @type nodevar: `str`
    @param ignore_initial: Do not include initial state info from TS.
        Enable this to mask absence of FTS initial states.
        Useful when initial states are specified in another way,
        e.g., directly augmenting the spec part.
    @type ignore_initial: `bool`
    @param receptive: if `True`, then add assumptions to
        ensure receptiveness at each node.

    @return: GR(1) formula representing `g`.
    @rtype: `GRSpec`
    """
    assert g.is_consistent()
    env_vars = {k: v for k, v in g.vars.iteritems()
                if k in g.env_vars}
    sys_vars = {k: v for k, v in g.vars.iteritems()
                if k not in g.env_vars}
    dom = _iter2var(g, nodevar)
    if g.owner == 'sys':
        sys_vars[nodevar] = dom
    elif g.owner == 'env':
        env_vars[nodevar] = dom
    else:
        raise ValueError('owner is "{owner}"'.format(owner=g.owner))
    dvars = dict(env_vars)
    dvars.update(sys_vars)
    # add primed copies -- same `dict` value
    p, _ = _prime_dict(dvars)
    dvars.update(p)
    evars = dict(env_vars)
    p, _ = _prime_dict(evars)
    evars.update(p)
    # convert to logic
    init = _init_from_ts(g.initial_nodes, nodevar, dvars, ignore_initial)
    tmp_init, nodepred = _node_var_trans(g, nodevar, dvars)
    if g.owner == 'sys':
        sys_init = init + tmp_init
        sys_safe = _sys_trans(g, nodevar, dvars)
        sys_safe += nodepred
        env_init = list()
        if receptive:
            env_safe = _env_trans_from_sys_ts(g, nodevar, dvars)
        else:
            env_safe = list()
    elif g.owner == 'env':
        sys_init = list()
        sys_safe = list()
        env_init = init + tmp_init
        env_safe = nodepred + _env_trans(g, nodevar, dvars)
    return GRSpec(
        sys_vars=sys_vars, env_vars=env_vars,
        env_init=env_init, sys_init=sys_init,
        env_safety=env_safe, sys_safety=sys_safe)


def _node_var_trans(g, nodevar, dvars):
    """Require variables to follow nodes according to labeling."""
    init = list()
    trans = list()
    # no AP labels ?
    if not dvars:
        return (init, trans)
    for u, d in g.nodes_iter(data=True):
        pre = _assign(nodevar, u, dvars)
        r = _to_action(d, dvars)
        if not r:
            continue
        # initial node vars
        init.append('!({pre}) || ({r})'.format(pre=pre, r=r))
        # transitions of node vars
        trans.append('(X (({pre}) -> ({r})))'.format(pre=pre, r=r))
    return (init, trans)


def _init_from_ts(initial_nodes, nodevar, dvars, ignore_initial=False):
    """Initial state, including enforcement of exactly one."""
    if ignore_initial:
        return list()
    if not initial_nodes:
        raise Exception(
            'FTS has no initial states.\n'
            'Enforcing this renders False the GR(1):\n'
            ' - guarantee if this is a system TS,\n'
            '   so the spec becomes trivially False.\n'
            ' - assumption if this is an environment TS,\n'
            '   so the spec becomes trivially True.')
    return [_disj(_assign(nodevar, u, dvars) for u in initial_nodes)]


def _sys_trans(g, nodevar, dvars):
    """Convert transition relation to GR(1) sys_safety."""
    logger.debug('modeling sys transitions in logic')
    sys_trans = list()
    for u in g.nodes_iter():
        pre = _assign(nodevar, u, dvars)
        # no successors ?
        if not g.succ.get(u):
            logger.debug('node: {u} is deadend !'.format(u=u))
            sys_trans.append('({pre}) -> (X False)'.format(pre=pre))
            continue
        post = list()
        for u, v, d in g.edges_iter(u, data=True):
            t = dict(d)
            t[_prime(nodevar)] = v
            r = _to_action(t, dvars)
            post.append(r)
        c = '({pre}) -> ({post})'.format(pre=pre, post=_disj(post))
        sys_trans.append(c)
    return sys_trans


def _env_trans_from_sys_ts(g, nodevar, dvars):
    """Convert environment actions to GR(1) env_safety.

    This constrains the actions available next to the environment
    based on the transition system.

    Purpose is to prevent env from blocking sys by purely
    picking a combination of actions for which sys has no outgoing
    transition from that state.
    """
    denv = {k: v for k, v in dvars.iteritems() if k in g.env_vars}
    env_trans = list()
    for u in g.nodes_iter():
        # no successor states ?
        if not g.succ.get(u):
            # nothing modeled for env, since sys has X(False) anyway
            # for action_type, codomain_map in env_action_ids.iteritems():
            # env_trans += [precond + ' -> X(' + s + ')']
            continue
        # collect possible next env actions
        c = set()
        for u, w, d in g.edges_iter(u, data=True):
            # TODO: syntactic over-approximation not applied,
            # so primed sys vars not filtered out here to
            # derive guards
            t = _to_action(d, denv)
            if not t:
                continue
            c.add(t)
        # no next env actions ?
        if not c:
            continue
        post = _disj(c)
        pre = _assign(nodevar, u, dvars)
        env_trans.append('(({pre}) -> ({post}))'.format(pre=pre, post=post))
    return env_trans


def _env_trans(g, nodevar, dvars):
    """Convert environment transitions to GR(1) safety assumption.

    @type g: `networkx.MutliDigraph`
    @param nodevar: name of variable representing current node
    @type nodevar: `str`
    @type dvars: `dict`
    """
    env_trans = list()
    for u in g.nodes_iter():
        pre = _assign(nodevar, u, dvars)
        # no successors ?
        if not g.succ.get(u):
            env_trans.append('{pre} -> X(False)'.format(pre=pre))
            warnings.warn(
                'Environment dead-end found.\n'
                'If sys can force env to dead-end,\n'
                'then GR(1) assumption becomes False,\n'
                'and spec trivially True.')
            continue
        post = list()
        sys = list()
        for u, v, d in g.out_edges_iter(u, data=True):
            # primed sys vars cannot appear in Mealy games
            # but can in Moore games.
            # Checked later, in `spec`.
            # Don't constrain the irrelevant.
            # action
            t = dict(d)
            t[_prime(nodevar)] = v
            r = _to_action(t, dvars)
            post.append(r)
            # what sys vars ?
            t = {k: v for k, v in d.iteritems()
                 if k not in g.env_vars}
            r = _to_action(t, dvars)
            sys.append(r)
        # avoid sys winning env by blocking all edges
        post.append(_conj_neg(sys))
        env_trans.append('({pre}) -> ({post})'.format(
            pre=pre, post=_disj(post)))
    return env_trans


def _to_action(d, dvars):
    """Return `str` conjoining assignments and `"formula"` in `d`.

    @param d: (partial) mapping from variables in `dvars`
        to values in their range, defined by `dvars`
    @type d: `dict`
    @type dvars: `dict`
    """
    c = list()
    if 'formula' in d:
        c.append(d['formula'])
    for k, v in d.iteritems():
        if k not in dvars:
            continue
        s = _assign(k, v, dvars)
        c.append(s)
    return _conj(c)


def _assign(k, v, dvars):
    """Return `str` of equality of variable `k` to value `v`.

    @type k: `str`
    @type v: `str` or `int`
    @type dvars: `dict`
    """
    dom = dvars[k]
    if isinstance(dom, tuple):
        s = '{k} = {v}'.format(k=k, v=v)
    elif isinstance(dom, (set, list)):
        s = '{k} = "{v}"'.format(k=k, v=v)
    elif dom in {'bool', 'boolean'}:
        s = '{k} <-> {v}'.format(k=k, v=v)
    else:
        raise Exception('domain is: {dom}'.format(dom=dom))
    return _pstr(s)


def build_dependent_var_table(fts, statevar):
    """Return a `dict` of substitution rules for dependent variables.

    The dependent variables in a transition system are the
    atomic propositions that are used to label states.

    They are "dependent" because their values are completely
    determined by knowledge of the current state in the
    transition system.

    The returned substitutions can be used

    @type fts: `TransitionSystem`

    @param statevar: name of variable used for the current state
        For example if it is 'loc', then the states
        `'s0', 's1'` are mapped to:

        ```
          {'s0': '(loc = "s0")',
           's1': '(loc = "s1")'}
        ```

    @type state_ids: `dict`

    @rtype: `{'p': '((loc = "s1") | (loc = "s2") | ...)', ...}`
        where:

          - `'p'` is a proposition in `fts.atomic_propositions`
          - the states "s1", "s2" are labeled with `'p'`
          - `loc` is the string variable used for the state of `fts`.
    """
    raise NotImplementedError('under maintenance')
    state_ids, __ = _iter2var(fts, variables=dict(), statevar=statevar,
                              bool_states=False, must='xor')
    ap2states = map_ap_to_states(fts)
    return {k: _disj(state_ids[x] for x in v)
            for k, v in ap2states.iteritems()}


def map_ap_to_states(fts):
    """For each proposition find the states labeled with it.

    @type fts: `TransitionSystem`

    @rtype: `{'p': s, ...}` where `'p'` a proposition and
        `s` a set of states in `fts`.
    """
    table = {p: set() for p in fts.vars}
    for u in fts:
        for p in fts.node[u]['ap']:
            table[p].add(u)
    return table


def synthesize_many(specs, ts=None, ignore_init=None, solver='gr1c'):
    """Synthesize from logic specs and multiple transition systems.

    The transition systems are composed synchronously, i.e.,
    they all have to take a transition at each time step.
    The synchronous composition is obtained by taking the
    conjunction of the formulas describing each transition system.

    The states of each transition system can be either:

      - all integers, or
      - all strings

    In either case the transition system state will be
    represented in logic with a single variable,
    that ranges over a finite set of integers or strings, respectively.

    The keys of `ts` are used to name each state variable.
    So the logic formula for `ts['name']` will be `'name'`.

    Who controls this state variable is determined from
    the attribute `TransitionSystem.owner` that can take the values:

      - `'env'`
      - `'sys'`

    For example:
    ```
    ts.add_nodes_from(xrange(4))
    ts['door'].owner = 'env'
    ```

    will result in a logic formula with
    an integer variable `'door'`
    controlled by the environment and
    taking values over `{0, 1, 2, 3}`.

    The example:
    ```
    ts.add_nodes_from(['a', 'b', 'c'])
    ts['door'].owner = 'sys'
    ```

    will instead result in a string variable `'door'`
    controlled by the system and taking
    values over `{'a', 'b', 'c'}`.

    @type specs: `GRSpec`

    @type ts: `dict` of `TransitionSystem`

    @type ignore_init: `set` of keys from `ts`

    @param solver: `'gr1c'` or `'jtlv'`
    @type solver: `str`
    """
    assert isinstance(ts, dict)
    for name, t in ts.iteritems():
        ignore = name in ignore_init
        statevar = name
        specs |= sys_to_spec(t, statevar, ignore)
    if solver == 'gr1c':
        ctrl = gr1c.synthesize(specs)
    elif solver == 'jtlv':
        ctrl = jtlv.synthesize(specs)
    else:
        raise Exception('Unknown solver: ' + str(solver) + '. '
                        'Available solvers: "jtlv" and "gr1c"')
    try:
        logger.debug('Mealy machine has: n = ' +
                     str(len(ctrl)) + ' states.')
    except:
        logger.debug('No Mealy machine returned.')
    # no controller found ?
    # counterstrategy not constructed by synthesize
    if not isinstance(ctrl, MealyMachine):
        return None
    remove_deadends(ctrl)
    return ctrl


def synthesize(option, specs, env=None, sys=None,
               ignore_env_init=False, ignore_sys_init=False,
               rm_deadends=True):
    """Function to call the appropriate synthesis tool on the specification.

    The states of the transition system can be either:

      - all integers, or
      - all strings

    For more details of how the transition system is represented in
    logic look at `synthesize_many`.


    Beware!
    =======
    This function provides a generic interface to a variety
    of routines.  Being under active development, the types of
    arguments supported and types of objects returned may change
    without notice.


    @param option: Magic string that declares what tool to invoke,
        what method to use, etc.  Currently recognized forms:

          - `"gr1c"`: use gr1c for GR(1) synthesis via `interfaces.gr1c`.
          - `"jtlv"`: use JTLV for GR(1) synthesis via `interfaces.jtlv`.
    @type specs: `spec.GRSpec`
    @param env: A transition system describing the environment:

            - states controlled by environment
            - input: `vars` not in `env_vars`
            - output: `env_vars`
            - initial states constrain the environment

        This constrains the transitions available to
        the environment, given the outputs from the system.
    @type env: `TransitionSystem`
    @param sys: A transition system describing the system:

            - states controlled by the system
            - input: `env_vars`
            - output: `vars` not in `env_vars`
            - initial states constrain the system
    @type sys: `TransitionSystem`
    @param ignore_sys_init: Ignore any initial state information
        contained in env.
    @type ignore_sys_init: `bool`
    @param ignore_env_init: Ignore any initial state information
        contained in sys.
    @type ignore_env_init: `bool`
    @param rm_deadends: return a strategy that contains no terminal states.
    @type rm_deadends: `bool`

    @return: If spec is realizable,
        then return a Mealy machine implementing the strategy.
        Otherwise return `None`.
    @rtype: `MealyMachine` or `None`
    """
    specs = _spec_plus_sys(specs, env, sys, ignore_env_init,
                           ignore_sys_init)
    if option == 'gr1c':
        strategy = gr1c.synthesize(specs)
    elif option == 'jtlv':
        strategy = jtlv.synthesize(specs)
    else:
        raise Exception('Undefined synthesis option. ' +
                        'Current options are "jtlv" and "gr1c"')
    ctrl = strategy2mealy(strategy, specs)
    try:
        logger.debug('Mealy machine has: n = ' +
                     str(len(ctrl)) + ' states.')
    except:
        logger.debug('No Mealy machine returned.')
    # no controller found ?
    # exploring unrealizability with counterexamples or other means
    # can be done by calling a dedicated other function, not this
    if not isinstance(ctrl, MealyMachine):
        return None
    if rm_deadends:
        remove_deadends(ctrl)
    return ctrl


def is_realizable(option, specs, env=None, sys=None,
                  ignore_env_init=False, ignore_sys_init=False):
    """Check realizability.

    For details, see `synthesize`.
    """
    specs = _spec_plus_sys(
        specs, env, sys,
        ignore_env_init, ignore_sys_init)
    if option == 'gr1c':
        r = gr1c.check_realizable(specs)
    elif option == 'jtlv':
        r = jtlv.check_realizable(specs)
    else:
        raise Exception('Undefined synthesis option. ' +
                        'Current options are "jtlv" and "gr1c"')
    if r:
        logger.debug('is realizable')
    else:
        logger.debug('is not realizable')
    return r


def _spec_plus_sys(specs, env, sys, ignore_env_init, ignore_sys_init):
    if sys is not None:
        assert sys.owner == 'sys'
        if hasattr(sys, 'state_varname'):
            statevar = sys.state_varname
        else:
            logger.info('sys.state_varname undefined. '
                        'Will use the default variable name: "loc".')
            statevar = 'loc'
        sys_formula = sys_to_spec(sys, statevar, ignore_sys_init)
        specs = specs | sys_formula
        logger.debug('sys TS:\n' + str(sys_formula.pretty()) + _hl)
    if env is not None:
        assert env.owner == 'env'
        if hasattr(env, 'state_varname'):
            statevar = sys.state_varname
        else:
            logger.info('env.state_varname undefined. '
                        'Will use the default variable name: "eloc".')
            statevar = 'eloc'
        env_formula = sys_to_spec(env, statevar, ignore_env_init)
        specs = specs | env_formula
        logger.debug('env TS:\n' + str(env_formula.pretty()) + _hl)
    logger.info('Overall Spec:\n' + str(specs.pretty()) + _hl)
    return specs


def strategy2mealy(A, spec):
    """Convert strategy to Mealy transducer.

    `A` is the contraction of the (deterministic) game graph
    that is the strategy.

    @param A: strategy
    @type A: `networkx.DiGraph`
    @type spec: `GRSpec`

    @return: `MealyMachine`
    """
    logger.info('converting strategy (compact) to Mealy machine')
    if not A:
        raise Exception(
            'Empty graph returned as synthesized strategy !\n'
            'Is your design perhaps trivially realizable ?\n'
            '(i.e., has false assumption ?)')
    env_vars = spec.env_vars
    sys_vars = spec.sys_vars
    mach = MealyMachine()
    inputs = machines.create_machine_ports(env_vars)
    mach.inputs.update(inputs)
    outputs = machines.create_machine_ports(sys_vars)
    mach.outputs.update(outputs)
    str_vars = {
        k: v for k, v in env_vars.iteritems()
        if isinstance(v, list)}
    str_vars.update({
        k: v for k, v in sys_vars.iteritems()
        if isinstance(v, list)})
    mach.add_nodes_from(A)
    # transitions labeled with I/O
    for u in A:
        for v in A.successors_iter(u):
            d = A.node[v]['state']
            d = _int2str(d, str_vars)
            mach.add_edge(u, v, **d)

            logger.info('node: {v}, state: {d}'.format(v=v, d=d))
    # special initial state, for first reaction
    initial_state = 'Sinit'
    mach.add_node(initial_state)
    mach.initial_nodes.add(initial_state)
    # fix an ordering for keys
    # because tuple(dict.iteritems()) is not safe:
    # https://docs.python.org/2/library/stdtypes.html#dict.items
    u = next(iter(A))
    keys = A.node[u]['state'].keys()
    # to store tuples of dict values for fast search
    isinit = spec.compile_init(no_str=True)
    # Mealy reaction to initial env input
    init_valuations = set()
    tmp = dict()
    for u, d in A.nodes_iter(data=True):
        var_values = d['state']
        vals = tuple(var_values[k] for k in keys)
        logger.debug(
            'machine vertex: {u}, has var values: {v}'.format(
                u=u, v=var_values))
        # already an initial valuation ?
        if vals in init_valuations:
            continue
        # add edge: Sinit -> u ?
        tmp.update(var_values)
        if not eval(isinit, tmp):
            continue
        label = _int2str(var_values, str_vars)
        mach.add_edge(initial_state, u, **label)
        # remember variable values to avoid
        # spurious non-determinism wrt the machine's memory
        init_valuations.add(vals)
        logger.debug('found initial state: {u}'.format(u=u))
    if mach.succ.get('Sinit'):
        return mach
    import pprint
    raise Exception(
        'The machine obtained from the strategy '
        'does not have any initial states !\n'
        'The strategy is:\n'
        'vertices:' + pprint.pformat(A.nodes(data=True)) + 2 * '\n' +
        'edges:\n' + str(A.edges()) + 2 * '\n' +
        'and the machine:\n' + str(mach) + 2 * '\n' +
        'and the specification is:\n' + str(spec.pretty()) + 2 * '\n')


def _int2str(label, str_vars):
    """Replace integers with string values for string variables.

    @param label: mapping from variable names, to integer (as strings)
    @type label: `dict`
    @param str_vars: mapping that defines those variables that
        should be converted from integer to string variables.
        Each variable is mapped to a list of strings that
        comprise its range. This list defines how integer values
        correspond to string literals for that variable.
    @type str_vars: `dict`

    @rtype: `dict`
    """
    label = dict(label)
    label.update({k: str_vars[k][int(v)]
                 for k, v in label.iteritems()
                 if k in str_vars})
    return label


def mask_outputs(machine):
    """Erase outputs from each edge where they are zero."""
    for u, v, d in machine.edges_iter(data=True):
        for k in d:
            if k in machine.outputs and d[k] == 0:
                d.pop(k)


def determinize_machine_init(mach, init_out_values=None):
    """Return a determinized copy of `mach` with given initial outputs.

    The transducers produced by synthesis can have multiple
    initial output valuations as possible reactions to a
    given input valuation.

    Possible reasons for this are:

      1. the system does not have full control over its initial state.
        For example the option `"ALL_INIT"` of `gr1c`.

      2. the strategy returned by the solver has multiple
        vertices that satisfy the initial conditions.

    Case 1
    ======
    Requires an initial condition to be specified for
    each run of the transducer, because the transducer does
    not have full freedom to pick the initial output values.

    Note that solver options like `"ALL_INIT"`
    assume that the system has no control initially.
    Any output valuation that satisfies the initial
    conditions can occur.

    However, this may be too restrictive.
    The system may have control over the initial values of
    some outputs, but not others.

    For the outputs it can initially control,
    the non-determinism resulting from synthesis is redundancy
    and can be removed arbitrarily, as in Case 2.

    Case 2
    ======
    The function `strategy2mealy` returns a transducer that
    for each initial input valuation,
    for each initial output valuation,
    reacts with a unique transition.

    But this can yield multile reactions to a single input,
    even for solver options like `"ALL_ENV_EXIST_SYS_INIT"` for `gr1c`.
    The reason is that there can be multiple strategy vertices
    that satisfy the initial conditions, but the solver
    included them not because they are needed as initial reactions,
    but to be visited later by the strategy.

    These redundant initial reactions can be removed,
    and because the system has full control over their values,
    they can be removed in an arbitrary manner,
    keeping only a single reaction, for each input valuation.

    Algorithm
    =========
    Returns a deterministic transducer.
    This means that at each transducer vertex,
    for each input valuation,
    there is only a single reaction (output valuation) available.

    The non-determinism is resolved for the initial reaction
    by ensuring the outputs given in `init_out_values`
    take those values.
    The remaining outputs are determinized arbitrarily.

    See also
    ========
    `synthesize`, `strategy2mealy`

    @param mach: possibly non-deterministic transducer,
        as produced, for example, by `synthesize`.
    @type mach: `MealyMachine`

    @param init_out_values: mapping from output ports that
        the system cannot control initially,
        to the initial values they take in this instance of the game.
    @type init_out_values: `dict`

    @rtype: `MealyMachine`
    """
    mach = copy.deepcopy(mach)
    if init_out_values is None:
        init_out_values = dict()
    '''determinize given outputs (uncontrolled)'''
    # restrict attention to given output ports only
    given_ports = tuple(k for k in mach.outputs if k in init_out_values)
    rm_edges = set()
    for i, j, key, d in mach.edges_iter(['Sinit'], data=True, keys=True):
        for k in given_ports:
            if d[k] != init_out_values[k]:
                rm_edges.add((i, j, key))
                break
    mach.remove_edges_from(rm_edges)
    '''determinize arbitrarily any remnant non-determinism'''
    # input valuations already seen
    # tuples of values used for efficiency (have __hash__)
    possible_inputs = set()
    # fix a key order
    inputs = tuple(k for k in mach.inputs)
    rm_edges = set()
    for i, j, key, d in mach.edges_iter(['Sinit'], data=True, keys=True):
        in_values = tuple(d[k] for k in inputs)
        # newly encountered input valuation ?
        if in_values not in possible_inputs:
            possible_inputs.add(in_values)
            continue
        else:
            rm_edges.add((i, j, key))
    mach.remove_edges_from(rm_edges)
    return mach
