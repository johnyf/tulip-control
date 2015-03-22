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
"""Algorithms on Kripke structures and Automata"""
import logging
logger = logging.getLogger(__name__)
from tulip.interfaces import ltl2ba as ltl2baint
from tulip.transys.automata import Automaton


_hl = 40 * '-'
# build parser once only
parser = ltl2baint.Parser()


def ltl2ba(formula):
    """Convert LTL formula to Buchi Automaton using ltl2ba.

    @type formula: `str(formula)` must be admissible ltl2ba input

    @return: Buchi automaton whose edges are annotated
        with Boolean formulas as `str`
    @rtype: [`Automaton`]
    """
    ltl2ba_out = ltl2baint.call_ltl2ba(str(formula))
    symbols, g, initial, accepting = parser.parse(ltl2ba_out)
    ba = Automaton('Buchi', alphabet=symbols)
    ba.add_nodes_from(g)
    ba.add_edges_from(g.edges_iter(data=True))
    ba.initial_nodes = initial
    ba.accepting_sets = accepting
    logger.info('Resulting automaton:\n\n{ba}\n'.format(ba=ba))
    return ba
