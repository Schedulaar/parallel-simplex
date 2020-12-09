import random


for n in [10, 100, 1000, 1000, 10000]:
    random.seed("Hello World!")
    fA = open(str(n) + "-A.csv", "w")
    fb = open(str(n) + "-b.csv", "w")
    fc = open(str(n) + "-c.csv", "w")
    for i in range(n):
        for j in range(n):
            fA.write('{:.4f}'.format(random.random()) + ",")
        fb.write('{:.4f}'.format(random.random()))
        fc.write('{:.4f}'.format(random.random()))
        if i < n-1:
            fA.write("\n")
            fb.write("\n")
            fc.write("\n")
    fA.close()
    fb.close()
    fc.close()
