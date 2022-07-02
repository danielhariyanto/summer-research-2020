import math
import numpy as np
import matplotlib.pyplot as plt

"""
FORMULAS:
βF = βFref + βFbond

βFbond = N*(math.log((1-pb)**f)+f/2*pb)
βFref = N*(math.log(density)-1+B2*density**2)
"""
"""
All βF and bond probability pb are dependent on temperature T and density ρ:
βF(T,ρ); βFref(T,ρ); βFbond(T,ρ); pb(T,ρ)
"""
"""
Fref(T,ρ) is the free energy of the reference state (the state in which no bonds are present)
Fbond(T,ρ) is the free-energy contribution due to the bonding between distinct tetramers
bond probability pb is the fraction of formed bonds over the total number of possible bonds, i.e., the probability that one arm of a tetramer is engaged in a bond

N is the number of particles
f is the particle valence
B2: We then ﬁt the density dependence of the osmotic pressure with the virial expression, estimating in this way the second virial coefﬁcient B2. We ﬁnd B2 = 2100±190nm6
"""

def main():
    N = 100
    pb = [0.3,0.4,0.5,0.6,0.7,0.8]
    f = [3,4,6,8,10,12]
    density = 20 #mg/ml
    B2 = 2100 #±190 nm**6

    #βFbond = calculateβFbond(N,pb,f)
    #βFref = calculateβFref(N,density,B2)
    #βF = calculateβF(βFbond,βFref)

    βFbondforpb = []
    for i in f:
        x = []
        for j in pb:
            y = calculateβFbond(N,j,i)
            x.append(y)
        βFbondforpb.append(x)

    calculateβFbondforpb(pb,f,βFbondforpb)

    βFbondforf = []
    for i in pb:
        x = []
        for j in f:
            y = calculateβFbond(N,i,j)
            x.append(y)
        βFbondforf.append(x)

    calculateβFbondforf(pb,f,βFbondforf)


def calculateβFbond(N,pb,f):
    βFbond = N*(math.log((1-pb)**f)+f/2*pb)
    return βFbond

def calculateβFref(N,density,B2):
    βFref = N*(math.log(density)-1+B2*density**2)
    return βFref

def calculateβF(βFbond,βFref):
    βF = βFref + βFbond
    return βF

def calculateβFbondforpb(pb,f,βFbondforpb):
    plt.plot(pb,βFbondforpb[0],label=str(f[0]),color='khaki',marker='o',markersize=1)
    plt.plot(pb,βFbondforpb[1],label=str(f[1]),color='gold',marker='o',markersize=1)
    plt.plot(pb,βFbondforpb[2],label=str(f[2]),color='goldenrod',marker='o',markersize=1)
    plt.plot(pb,βFbondforpb[3],label=str(f[3]),color='darkgoldenrod',marker='o',markersize=1)
    plt.plot(pb,βFbondforpb[4],label=str(f[4]),color='saddlebrown',marker='o',markersize=1)
    plt.plot(pb,βFbondforpb[5],label=str(f[5]),color='brown',marker='o',markersize=1)

    plt.xlabel("bond probability")
    plt.ylabel("Helmholtz free energy")
    plt.title("Estimated Helmholtz free energy values for various # of arms and bond probabilities")
    plt.legend(title='# of arms')

    plt.show()

def calculateβFbondforf(pb,f,βFbondforf):
    plt.plot(f,βFbondforf[0],label=str(pb[0]),color='khaki',marker='o',markersize=1)
    plt.plot(f,βFbondforf[1],label=str(pb[1]),color='gold',marker='o',markersize=1)
    plt.plot(f,βFbondforf[2],label=str(pb[2]),color='goldenrod',marker='o',markersize=1)
    plt.plot(f,βFbondforf[3],label=str(pb[3]),color='darkgoldenrod',marker='o',markersize=1)
    plt.plot(f,βFbondforf[4],label=str(pb[4]),color='saddlebrown',marker='o',markersize=1)
    plt.plot(f,βFbondforf[5],label=str(pb[5]),color='brown',marker='o',markersize=1)

    plt.xlabel("# of arms")
    plt.ylabel("Helmholtz free energy")
    plt.title("Estimated Helmholtz free energy values for various # of arms and bond probabilities")
    plt.legend(title='bond probabilities')

    plt.show()

main()
