import sympy as sp

# Symbolic expansion of the unperturbed diagonal selector for the s=2 kernel.
t, x, y = sp.symbols('t x y', real=True)

def state_sum(sig, tau):
    total = 0
    for eps in (1, -1):
        for eta in (1, -1):
            vx = eta - eps
            vy = eta*tau - eps*sig
            z = 2*t*(x*vx + y*vy) + t**2*(vx**2 + vy**2)
            total += sp.series(1/(1+z), t, 0, 7).removeO()
    return sp.expand(total)

F = {(s,u): state_sum(s,u) for s in (1,-1) for u in (1,-1)}
J = sp.expand((F[(1,1)] - F[(1,-1)] - F[(-1,1)] + F[(-1,-1)])/4)
h = sp.expand((F[(1,1)] + F[(1,-1)] - F[(-1,1)] - F[(-1,-1)])/4)

J4 = sp.expand(J.coeff(t,4)).subs(x**2 + y**2, 1)
h2 = sp.expand(h.coeff(t,2))
h4 = sp.expand(h.coeff(t,4)).subs(x**2 + y**2, 1)

print('J t^4 coefficient before x^2+y^2 simplification:', sp.expand(J.coeff(t,4)))
print('Expected angular form: 1536*x^2*y^2 - 160 = 384*sin(2theta)^2 - 160')
print('h t^2 coefficient:', h2, '= 16*sin(2theta)')
print('h t^4 coefficient before simplification:', sp.expand(h.coeff(t,4)))
print('Expected simplified h t^4: 256*x*y = 128*sin(2theta)')

z5 = float(sp.zeta(5)-1)
q6_bound = float(sp.zeta(3)**2/8)
p6 = float(sp.zeta(5)-sp.zeta(6))
local = 225*(0.037 + 2*0.181 + 0.020)
obj_halo = 225*z5
print('zeta(5)-1 =', z5)
print('Q6 upper bound zeta(3)^2/8 =', q6_bound)
print('P6 = zeta(5)-zeta(6) =', p6)
print('conservative local tree budget =', local)
print('objective halo bound =', obj_halo)
print('distinguished leading objective coefficient =', 224/8)
print('conservative objective lower bound =', 223/8 - obj_halo)
