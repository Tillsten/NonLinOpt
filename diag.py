import sympy
#from itertools import accumulate
from sympy.parsing.sympy_parser import parse_expr

greek_letters = {'a': r'\alpha', 'b': r'\beta',
                 'c': r'\gamma', 'd': r'\delta',
                 'e': r'\epsilon'}
for i in range(10):
    greek_letters[str(i)] = str(i)

import matplotlib.pyplot as plt

def m(s):
    """Puts $ around a string"""
    return r'$' + s + r'$'


def draw_states(pos, side, trans):
    """Draws the states"""
    off = -0.2 if side else 0.2
    plt.text(side + off, pos + 1, m(greek_letters[trans[-1]]),
        fontsize=20, horizontalalignment='center')


def draw_arrows(pos, side, direction, om):
    dx, dy = 0.5, 1
    npos = pos + 1
    if isinstance(om, int):
        txt = m(r'\omega_' + str(om))
    else:
        txt = m(r'k_{' + str(om) + r'}')
    if (side == 1 and direction == 'in') or (side == 0 and direction == 'out'):
        txt = r'$-' + txt[1:]
    else:
        txt = r'$+' + txt[1:]
    if direction == 'in':
        if side == 1:
            dx *= -1
        plt.arrow(side - dx, npos - dy, dx, dy,
            width=0.04, length_includes_head=True, lw=0)
        plt.text(side - dx, npos - dy / 2., txt, fontsize=20,
            horizontalalignment='center')
    else:
        if side == 1:
            dx *= -1
        plt.arrow(side, npos, -dx, dy, width=0.04, length_includes_head=True,
            lw=0)
        plt.text(side - dx, npos + dy / 3., txt, fontsize=20,
            horizontalalignment='center')


def mult_str(l):
    l = map(lambda x: x + '*', l)
    return l


class FeyDiag(object):
    """Feynman diagramm of non-linear optics"""

    def __init__(self, start_state, interactions):
        """
         A Feynman diagram is completely defined by the start state
         and the interactions.
        """
        self.start_state = start_state
        #self.om_to_trans = om_to_trans
        self.states = set(
            [c for ome, trans, di, side in interactions for c in trans])
        syms = sympy.symbols(self.states)
        #for i in syms:
        #    print i
        self.interactions = interactions
        self.left = [i for i in self.interactions if i[-1] == 0]
        self.right = [i for i in self.interactions if i[-1] == 1]


    def count_right(self):
        num = -1 ** (len(
            [' ' for om, trans, direction, side in self.interactions if
             side == 1]))
        if num == -1:
            return '-'
        else:
            return ''

    def formula(self):
        """
        Prints out the formula resulting of the interactions.
        All in the impulsive limit.
        """
        initial_pop = r'\rho_{%s%s}' % (
        greek_letters[self.start_state], greek_letters[self.start_state])
        left_tdms = [r'\mu_{%s}' % i[1] for i in self.left]
        right_tdms = [r'\mu_{%s}' % i[1] for i in self.right]
        return r'${0:>s}{1:>s}{2:>s}{3:>s}{4:>s}$'.format(
            str(self.count_right()), ''.join(left_tdms), initial_pop,
            ''.join(right_tdms), ''.join(self.propagators()))

    def formula_notex(self):
        from sympy.parsing.sympy_parser import parse_expr

        propagators = []
        for  i, (om, trans, direction, side) in enumerate(self.interactions):
            s = r'exp(-i*omega_%s*tau_%s)' % (trans, str(i))
            propagators.append(s)
        initial_pop = [r'p_%s%s' % (self.start_state, self.start_state)]
        left_tdms = [r'mu_%s' % i[1] for i in self.left]
        right_tdms = [r'mu_%s' % i[1] for i in self.right]
        f = initial_pop + left_tdms + right_tdms + propagators[:-1]
        return f

    def propagators(self):
        l = []
        for  i, (om, trans, direction, side) in enumerate(self.interactions):
            s = r'\exp(-i\omega_{%s}\tau_%s)' % (trans, str(i))
            l.append(s)
        return l

    def draw(self):
        height = len(self.interactions)
        plt.vlines(0, 0, height + 1, lw=5)
        plt.vlines(1, 0, height + 1, lw=5)
        plt.xlim(-1, 2)
        plt.ylim(-1, height + 1)
        #plt.figsize(*plt.figaspect(1))
        initial_pop = r'\rho_{%s%s}' % (self.start_state, self.start_state)
        plt.text(1 / 2., 0, m(initial_pop), fontsize=30,
            horizontalalignment='center')
        plt.text(1 / 2., -0.75, self.formula(), fontsize=15,
            horizontalalignment='center')
        plt.axis('off')
        pos = 0
        for  om, trans, direction, side in self.interactions:
            draw_arrows(pos, side, direction, om)
            draw_states(pos, side, trans)
            pos += 1


def sub_omega(formula):
    terms = formula.free_symbols
    for t in terms:
        if t.name.startswith('omega_'):
            a, b = t.name[-2], t.name[-1]
            ea = 'epsilon_' + a
            eb = 'epsilon_' + b
            formula = formula.subs(t,
                '(%s - %s)/hbar - i/Gamma_%s' % (eb, ea, a + b))
    return formula


def make_sympy(f_str):
    f_str = mult_str(f_str)
    return parse_expr(''.join(f_str)[:-1]).simplify()


def make_field(response):
    return response / sympy.Symbol('hbar')


def generate_nth_order(n):
    sol = []

    def visit(prev, actions, bra, ket):
        if bra - ket - 1 > actions:
            return
        a = actions - 1
        if actions > 0:
            visit(prev[:] + [
                (str(n - actions), '%s%s' % (bra, bra + 1), 'in', 0)], a,
                bra + 1, ket)
            visit(prev[:] + [
                (str(n - actions), '%s%s' % (ket, ket + 1), 'in', 1)], a, bra,
                ket + 1)
            if bra > 0:
                visit(prev[:] + [
                    (str(n - actions), '%s%s' % (bra, bra - 1), 'out', 0)], a,
                    bra - 1, ket)
            if ket > 0:
                visit(prev[:] + [
                    (str(n - actions), '%s%s' % (ket, ket - 1), 'out', 1)], a,
                    bra, ket - 1)
        elif actions == 0:
            if bra - ket == 1:
                sol.append(prev + [('sig', '%s%s' % (bra, bra - 1), 'out', 0)])

    visit([], n, 0, 0)
    return sol

print generate_nth_order(3)


def test():
    om = sympy.Symbol('omega_ab')
    ans = parse_expr('(-epsilon_a+epsilon_b)/hbar-i/Gamma_ab')
    first_order = FeyDiag(start_state=r'a',
        interactions=[('pu', 'ab', 'in', 0), ('pu', 'ab', 'out', 0)])
    f = first_order.formula_notex()
    f = make_sympy(f)
    f = sub_omega(f)
    e = make_field(f)
    det = sympy.integrate(e, ('tau_0', -sympy.oo, sympy.oo))

    sympy.pprint(det)

    assert(ans == sub_omega(om) )

#test()

