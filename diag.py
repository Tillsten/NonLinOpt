import sympy
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

def draw_arrows(pos, side, direction, om, last=False):
    ARROW_WIDTH = 0.04
    FONTSIZE = 20
    alpha = 0.3 if last else None
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
            width=ARROW_WIDTH, length_includes_head=True, lw=0)
        plt.text(side - dx, npos - dy / 2., txt, fontsize=FONTSIZE,
            horizontalalignment='center')
    else:
        if side == 1:
            dx *= -1
        plt.arrow(side, npos, -dx, dy, width=ARROW_WIDTH,
            length_includes_head=True,
            lw=0, alpha=alpha)
        plt.text(side - dx, npos + dy / 3., txt, fontsize=FONTSIZE,
            horizontalalignment='center')

def mult_str(l):
    l = map(lambda x: x + '*', l)
    return l

def sub_omega(formula):
    terms = formula.free_symbols
    for t in terms:
        if t.name.startswith('omega_'):
            a, b = t.name[-2], t.name[-1]
            ea = 'epsilon_' + a
            eb = 'epsilon_' + b
            formula = formula.subs(t,
                'omega_%s - I/Gamma_%s' % (a + b, a + b))
    return formula

def eq_mu(formula):
    terms = formula.free_symbols
    for t in terms:
        if t.name.startswith('mu_'):
            a, b = t.name[-2], t.name[-1]
            formula = formula.subs(t, 'mu_' + b + a)
    return formula

def make_sympy(f_str):
    f_str = mult_str(f_str)
    return parse_expr(''.join(f_str)[:-1]).simplify()

def make_field(response):
    return response / sympy.Symbol('hbar')

def generate_nth_order(n,start=(0,0)):
    sol = []
    def visit(prev, actions, bra=0, ket=0):
        if bra - ket - 1 > actions:
            return
        a = actions - 1        
        if actions > 0:
            visit(prev[:] + [
                (str(n - actions), (0,1), 0)], a,
                bra + 1, ket)
            visit(prev[:] + [
                (str(n - actions), (0,1), 1)], a, bra,
                ket + 1)
            if bra > 0:
                visit(prev[:] + [
                    (str(n - actions), (0,-1), 0)], a,
                    bra - 1, ket)
            if ket > 0:
                visit(prev[:] + [
                    (str(n - actions), (0,-1), 1)], a,
                    bra, ket - 1)
        elif actions == 0:
            if bra - ket == 1:
                sol.append([((0,),(0,))]+prev + [('sig', (0,-1), 0)])
            
    visit([], n, *start)
    return sol
    

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
    
#New format, example for a double oscilator:
bra, ket = (0,0),(0,0)
start_state = (bra, ket)
#format is name, delta bra(ket), ket or bar                       
example_interaction = ('pr',(1,1), 0) 
example_diagramm = [start_state,
                    (('pu'),(1,1),0), (('pu'),(1,-1), 0), 
                    (('pr'),(1,1),0), (('pr'),(1,-1),0)]

test_diag = [((0,),(0,)), (('pu'),(0,1),0), (('pu'),(0,1), 1), 
             (('pr'),(0,1),0), (('pr'),(0,-1),0)]

def draw(diag):
    """Draws the diagram"""
    height = len(diag)-1
    plt.vlines(0, 0, height + 1, lw=5)
    plt.vlines(1, 0, height + 1, lw=5)
    plt.xlim(-1, 2)
    plt.ylim(-1, height + 1)
    draw_state(-1, diag[0])
    f = m(sympy.printing.latex(formula_sympy(diag)))
    plt.text(1 / 2., -0.75, f, fontsize=25,
        horizontalalignment='center')    
    plt.axis('off')
    pos = 0
    state = diag[0]   
    for name, delta, side in diag[1:]:
        last = pos == len(diag) - 2
        state=apply_delta(state, delta, side)
        if sum(delta) > 0: 
            direction = 'in'
        else: 
            direction ='out'            
        draw_arrows(pos, side, direction, name, last)
        draw_state(pos, state)                
        pos += 1

def draw_state(pos, state):    
    s1, s2 = str(state[0]).strip('(),[]'), str(state[1]).strip('(),[]')
    plt.text(0.5, pos + 1, m(r'\left|  %s\rangle\langle%s\right|'%(s1,s2)),
        fontsize=25, horizontalalignment='center', verticalalignment='center')

def apply_delta(state, delta, side):
    t = list(state[side])
    idx, value = delta
    t[idx] += value
    if side==0:
        return tuple(t), state[1]
    else:
        return state[0], tuple(t)


def get_trans(state, delta, side):    
    old = state[side][delta[0]]
    new = old+delta[1]
    return old, new
    

def count_right(diag):
    return sum([i[-1] for i in diag[1:]])

def formula_sympy(diag):
    state = diag[0]
    propagators = []
    tdms = []
    for  i, (name, delta, side) in enumerate(diag[1:]):        
        om = get_trans(state, delta, side)
        state = apply_delta(state, delta, side)
        s = r'exp(-I*omega_%s%s*tau_%s)' % (om[0],om[1], str(i))
        propagators.append(s)           
        tdms.append('mu_%s%s'%(om[0],om[1]))        
    initial_pop = [r'p_00']
    
    f = initial_pop + tdms + propagators[:-1]
    f = make_sympy(f)
    f = eq_mu(f)
    f = sub_omega(f).simplify()    
    f = f.subs('tau_0',0)
    f = f.subs('tau_1',0)
    return f*(-1)**count_right(diag)

#draw(test_diag)
#print generate_nth_order(3)[2]
#draw(generate_nth_order(3)[2])

l=generate_nth_order(3)
for i in l:
    draw(i)
    figure()
        #pri
#test()

