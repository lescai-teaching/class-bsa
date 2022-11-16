import math
import matplotlib.pyplot as plt

def probCov(limit, coverage):
    prob = ( pow(coverage,limit) * math.exp(-coverage) ) / (math.factorial(limit))
    return(prob)


def sumProb(threshold, coverage):
    probtotal = 0
    for limit in range(0,threshold+1):
        probtotal = probtotal + probCov(limit, coverage)
    return(probtotal)


x = list(range(1,30))
y = list()
for cov in range(1,30):
    y.append(sumProb(5,cov))
plt.plot(x,y)
plt.show()