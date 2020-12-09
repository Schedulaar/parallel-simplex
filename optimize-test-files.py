from scipy.optimize import linprog
from numpy import genfromtxt
import time

print("Please enter matrix size nxn.")
n = int(input())

A = genfromtxt(str(n) + '-A.csv', delimiter=',')[:n, :n]
b = genfromtxt(str(n) + '-b.csv', delimiter='\n')
c = genfromtxt(str(n) + '-c.csv', delimiter='\n')

print("Starting...")
start = time.time()

res = linprog(-c, A, b, bounds=(0, None), method='revised simplex')

end = time.time()

print(res)
print("Finished with objective value " + str(res.fun))
print("This took " + str(end - start) + " seconds.")