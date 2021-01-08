from pysmps import smps_loader as loader
import sys

def truncate(f):
    """Truncates/pads a float f to n decimal places without rounding"""
    s = '{:f}'.format(f)
    s.strip("0")
    i, p, d = s.partition('.')
    return '.'.join([i, d.rstrip("0")])

result = loader.load_mps("netlib/" + sys.argv[1] + ".SIF")

(name, objective_name, row_names, col_names, _, types, c, A, rhs_names, rhs, bnd_names, bnd) = result

if len(rhs_names) != 1:
    raise Exception("Can't handle more than one b vector")

b = rhs[rhs_names[0]]

fA = open("netlib/" + name + "-A.csv", "w")
fb = open("netlib/" + name + "-b.csv", "w")
fc = open("netlib/" + name + "-c.csv", "w")
fshape = open("netlib/" + name + "-shape.csv", "w")

(m, n) = A.shape
actualM = m

for i in range(m):
    if types[i] == "L":
        for j in range(n):
            fA.write(truncate(A[i, j]) + ",")
        fb.write(truncate(b[i]))
        if i < m - 1:
            fA.write("\n")
            fb.write("\n")
    if types[i] == "G":
        for j in range(n):
            fA.write(truncate(-A[i, j]) + ",")
        fb.write(truncate(-b[i]))
        if i < m - 1:
            fA.write("\n")
            fb.write("\n")
    if types[i] == "E":
        for j in range(n):
            fA.write(truncate(A[i, j]) + ",")
        fb.write(truncate(b[i]))
        fA.write("\n")
        fb.write("\n")
        for j in range(n):
            fA.write(truncate(-A[i, j]) + ",")
        fb.write(truncate(-b[i]))
        if i < m - 1:
            fA.write("\n")
            fb.write("\n")
        actualM += 1

fshape.write('{}\n{}'.format(actualM, n))
fshape.close()

for j in range(n):
    fc.write(truncate(-c[j]))
    if j < n - 1:
        fc.write("\n")

fA.close()
fb.close()
fc.close()
