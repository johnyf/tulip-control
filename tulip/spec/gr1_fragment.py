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
#
"""Test if given formula belongs to an LTL fragment that
is convertible to deterministic Buchi Automata
(readily expressible in GR(1) ).

reference
=========
1. Andreas Morgenstern and Klaus Schneider,
   A LTL Fragment for GR(1)-Synthesis,
   in Proceedings First International Workshop on
   Interactions, Games and Protocols (iWIGP),
   Electronic Proceedings in Theoretical Computer Science (EPTCS),
   50, pp. 33--45, 2011,
   http://doi.org/10.4204/EPTCS.50.3
"""
from __future__ import absolute_import
import logging
logger = logging.getLogger(__name__)
from tulip import transys as trs
from tulip.spec import lexyacc, GRSpec
from tulip.spec import transformation as tx
from tulip.spec import ast as sast


def check(formula):
    """Parse formula string and create abstract syntax tree (AST)."""
    ast = lexyacc.parse(formula)
    nodes = {'gf', 'fg', 'g', 'f'}
    dfa.add_nodes_from(nodes)
    dfa.initial_nodes.add('gf')

    dfa = trs.automata.FiniteWordAutomaton(atomic_proposition_based=False,
                                           deterministic=True)

    dfa.alphabet |= {'!', 'W', 'U', 'G', 'F',
                     'U_left', 'U_right',
                     'W_left', 'W_right'}

    dfa.add_edge('gf', 'fg', letter='!')
    dfa.add_edge('fg', 'gf', letter='!')
    dfa.add_edge('g', 'f', letter='!')
    dfa.add_edge('f', 'g', letter='!')

    dfa.add_edge('gf', 'gf', letter='W')
    dfa.add_edge('gf', 'gf', letter='U_left')
    dfa.add_edge('gf', 'gf', letter='G')

    dfa.add_edge('fg', 'fg', letter='U')
    dfa.add_edge('fg', 'fg', letter='F')
    dfa.add_edge('fg', 'fg', letter='W_right')

    dfa.add_edge('gf', 'f', letter='U_right')
    dfa.add_edge('gf', 'f', letter='F')

    dfa.add_edge('fg', 'g', letter='W_left')
    dfa.add_edge('fg', 'g', letter='G')

    dfa.add_edge('g', 'g', letter='W')
    dfa.add_edge('g', 'g', letter='G')

    dfa.add_edge('f', 'f', letter='U')
    dfa.add_edge('f', 'f', letter='F')

    # plot tree automaton
    # dfa.save('dfa.pdf')

    # plot parse tree
    sast.dump_dot(ast, 'ast.dot')

    # sync product of AST with DFA,
    # to check acceptance
    Q = [(ast, 'gf')]
    while Q:
        s, q = Q.pop()
        logger.info('visiting: ' + str(s) + ', ' + str(q))

        if isinstance(s, sast.Unary):
            op = s.operator

            if op in {'!', 'G', 'F'}:
                t = _filter_edges(dfa, q, dict(letter=op))

                if not t:
                    raise Exception('not in fragment')

                qi, qj, w = t[0]

                Q.append((s.operand, qj))
            else:
                # ignore
                Q.append((s.operand, q))
        elif isinstance(s, sast.Binary):
            op = s.operator

            if op in {'W', 'U'}:
                t = _filter_edges(dfa, q, dict(letter=op))
                if t:
                    qi, qj, w = t[0]
                    Q.append((s.op_l, qj))
                    Q.append((s.op_r, qj))
                else:
                    t = _filter_edges(dfa, q, dict(letter=op + '_left'))

                    if not t:
                        raise Exception('not in fragment')

                    qi, qj, w = t[0]

                    Q.append((s.op_l, qj))

                    t = _filter_edges(dfa, q, dict(letter=op + '_right'))

                    if not t:
                        raise Exception('not in fragment')

                    qi, qj, w = t[0]

                    Q.append((s.op_r, qj))
            else:
                # ignore
                Q.append((s.op_l, q))
                Q.append((s.op_r, q))
        elif isinstance(s, sast.Var):
            print('reached var')

    return ast


def _filter_edges(g, u, attr):
    return [e for e, d in g.edges(u, data=True) if d == attr]


def stability_to_gr1(p, aux='aux'):
    """Convert C{<>[] p} to GR(1).

    Warning: This conversion is sound, but not complete.
    See p.2, U{[E10]
    <http://tulip-control.sourceforge.net/doc/bibliography.html#e10>}

    GR(1) form::

        !(aux) &&
        [](aux -> X aux) &&
        []<>(aux) &&

        [](aux -> p)

    @type p: str

    @param aux: name to use for auxiliary variable
    @type aux: str

    @rtype: L{GRSpec}
    """
    logging.warning(
        'Conversion of stability (<>[]p) to GR(1)' +
        'is sound, but NOT complete.'
    )

    a = aux
    a0 = a

    p = _paren(p)
    a = _paren(a)

    v = tx.check_var_name_conflict(p, a0)

    sys_vars = v | {a0}
    sys_init = {'!' + a}
    sys_safe = {a + ' -> ' + p,
                a + ' -> X ' + a}
    sys_prog = {a}

    return GRSpec(sys_vars=sys_vars, sys_init=sys_init,
                  sys_safety=sys_safe, sys_prog=sys_prog)


def response_to_gr1(p, q, aux='aux'):
    """Convert C{[](p -> <> q)} to GR(1).

    GR(1) form::

        []<>(aux) &&

        []( (p && !q) -> X ! aux) &&
        []( (! aux && !q) -> X ! aux)

    @type p: str

    @type q: str

    @param aux: name to use for auxiliary variable
    @type aux: str

    @rtype: L{GRSpec}
    """
    a = aux
    a0 = a

    p = _paren(p)
    q = _paren(q)
    a = _paren(a)

    s = p + ' -> <> ' + q
    v = tx.check_var_name_conflict(s, a0)

    sys_vars = v | {a0}
    # sys_init = {a}
    sys_safe = {
        '(' + p + ' && !' + q + ') -> X !' + a,
        '(!' + a + ' && !' + q + ') -> X !' + a
    }
    sys_prog = {a}

    return GRSpec(sys_vars=sys_vars,  # sys_init=sys_init,
                  sys_safety=sys_safe, sys_prog=sys_prog)


def eventually_to_gr1(p, aux='aux'):
    """Convert C{<> p} to GR(1).

    GR(1) form::

        !(aux) &&
        [](aux -> X aux) &&
        []<>(aux) &&

        []( (!p && !aux) -> X!(aux) )

    @type p: str

    @param aux: name to use for auxiliary variable
    @type aux: str

    @rtype: L{GRSpec}
    """
    a = aux
    a0 = a

    p = _paren(p)
    a = _paren(a)

    v = tx.check_var_name_conflict(p, a0)

    sys_vars = v | {a0}
    sys_init = {'!(' + a + ')'}
    sys_safe = {
        '(!' + p + ' && !' + a + ') -> X !' + a,
        a + ' -> X ' + a
    }
    sys_prog = {a}

    return GRSpec(sys_vars=sys_vars, sys_init=sys_init,
                  sys_safety=sys_safe, sys_prog=sys_prog)


def until_to_gr1(p, q, aux='aux'):
    """Convert C{p U q} to GR(1).

    GR(1) form::

        (!q -> !aux) &&
        [](q -> aux)
        [](aux -> X aux) &&
        []<>(aux) &&

        []( (!aux && X(!q) ) -> X!(aux) ) &&
        [](!aux -> p)

    @type p: str

    @param aux: name to use for auxiliary variable
    @type aux: str

    @rtype: L{GRSpec}
    """
    a = aux
    a0 = a

    p = _paren(p)
    q = _paren(q)
    a = _paren(a)

    s = p + ' && ' + q
    v = tx.check_var_name_conflict(s, a0)

    sys_vars = v | {a0}
    sys_init = {'!' + q + ' -> !' + a}
    sys_safe = {
        q + ' -> ' + a,
        '( (X !' + q + ') && !' + a + ') -> X !' + a,
        a + ' -> X ' + a,
        '(!' + a + ') -> ' + p
    }
    sys_prog = {a}

    return GRSpec(sys_vars=sys_vars, sys_init=sys_init,
                  sys_safety=sys_safe, sys_prog=sys_prog)


def _paren(x):
    return '({x})'.format(x=x)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    s = '(a U b) && []a && <>a && <>a && []<>(<>z)'
    parsed_formula = check(s)

    print('Parsing result: ' + str(parsed_formula))
