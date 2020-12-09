import random

random.seed("Hello World!")

for n in [10, 100, 1000, 1000, 10000, 100000, 1000000, 10000000, 100000000]:
    fA = open(str(n) + "-A.csv", "w")
    fb = open(str(n) + "-b.csv", "w")
    fc = open(str(n) + "-c.csv", "w")
    for i in range(n):
        for j in range(n):
            fA.write('{:.4f}'.format(random.random()) + ",")
        fA.write("\n")
        fb.write('{:.4f}'.format(random.random()) + "\n")
        fc.write('{:.4f}'.format(random.random()) + "\n")
    fA.close()
    fb.close()
    fc.close()
