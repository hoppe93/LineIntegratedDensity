#!/usr/bin/env python3

import sympy

xi, xi1, zi, zi1 = sympy.symbols('xi xi1 zi zi1')
t, l = sympy.symbols('t l')
x0, y0, z0 = sympy.symbols('x0 y0 z0')
nx, ny, nz = sympy.symbols('nx ny nz')

eq1 = zi + t*(zi1 - zi) - z0 - l*nz
sol1 = sympy.solve(eq1, t)

#print(sol1)

eq2 = (xi+((z0-zi+l*nz)/(zi1-zi))*(xi1-xi))**2 - (x0+l*nx)**2 - (y0+l*ny)**2
#sol2 = sympy.roots(eq2, l)

poly = sympy.poly(sympy.expand(eq2), l)
coeffs = poly.coeffs()

a0 = coeffs[2]
a1 = coeffs[1]
a2 = coeffs[0]

print(f'{len(coeffs)} coefficients')
print(':: a2')
print(f'  {a2}\n')
print(':: a1')
print(f'  {a1}\n')
print(':: a0')
print(f'  {a0}\n')

my_a0 = xi**2 - x0**2 - y0**2 + 2*xi*(z0-zi)*(xi1-xi)/(zi1-zi) + (z0-zi)**2*((xi1-xi)/(zi1-zi))**2
my_a1 = 2*xi*nz*(xi1-xi)/(zi1-zi) + 2*nz*(z0-zi)*((xi1-xi)/(zi1-zi))**2 - 2*(x0*nx + y0*ny)
my_a2 = nz**2*((xi1-xi)/(zi1-zi))**2 - nx**2 - ny**2

print(f'a0 = my_a0?  {sympy.simplify(my_a0-a0)==0}')
print(f'a1 = my_a1?  {sympy.simplify(my_a1-a1)==0}')
print(f'a2 = my_a2?  {sympy.simplify(my_a2-a2)==0}')

