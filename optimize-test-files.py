from scipy.optimize import linprog
from numpy import genfromtxt

n = 1000

A = genfromtxt(str(n) + '-A.csv', delimiter=',')[:n, :n]
b = genfromtxt(str(n) + '-b.csv', delimiter='\n')
c = genfromtxt(str(n) + '-c.csv', delimiter='\n')

res = linprog(-c, A, b, bounds=(0, None))
print(res)