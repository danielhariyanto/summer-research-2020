import math
import numpy as np
import matplotlib.pyplot as plt

def main():
    sequence = ['GC','GCGC','GCTAGC','GCTCGAGC','GCTCAATTGAGC','GCTCAACGGCTTGAGC']

    temperature = np.arange(20,80,0.1).tolist()

    Arms = [24,12,8,6,4,3]

    GValues = []
    for s in sequence:
        GValue = []
        for t in temperature:
            GValue.append(calculateG(s,t))
        GValues.append(GValue)

    HValues = [calculateH(sequence[0]),calculateH(sequence[1]),calculateH(sequence[2]),calculateH(sequence[3]),calculateH(sequence[4]),calculateH(sequence[5])]
    SValues = [calculateS(sequence[0]),calculateS(sequence[1]),calculateS(sequence[2]),calculateS(sequence[3]),calculateS(sequence[4]),calculateS(sequence[5])]

    GValuesStickyNucleotides = []
    for g in range(len(GValues)):
        GValueStickyNucleotides = []
        for h in GValues[g]:
            GValueStickyNucleotides.append(h*Arms[g])
        GValuesStickyNucleotides.append(GValueStickyNucleotides)

    HValuesStickyNucleotides = []
    for h in range(len(HValues)):
        HValuesStickyNucleotides.append(HValues[h]*Arms[h])

    SValuesStickyNucleotides = []
    for s in range(len(SValues)):
        SValuesStickyNucleotides.append(SValues[s]*Arms[s])


    """----------------------------ΔG----------------------------"""
    """Estimated ΔG values for corresponding transition temperatures"""
    TvsG(GValues,temperature)

    """Estimated transition temperatures for corresponding ΔG values"""
    #GvsT(GValues,temperature)

    """Estimated ΔG values of sticky nucleotides for corresponding transition temperatures"""
    TvsGwithStickyNucleotides(GValuesStickyNucleotides,temperature)

    """----------------------------ΔH----------------------------"""
    """Estimated ΔH values for corresponding # of sticky ends"""
    StickyEndvsH(HValues,sequence)
    HwithStickyNucleotides(HValuesStickyNucleotides,sequence)

    """----------------------------ΔS----------------------------"""
    """Estimated ΔS values for corresponding # of sticky ends"""
    StickyEndvsS(SValues,sequence)
    SwithStickyNucleotides(SValuesStickyNucleotides,sequence)


"""Calculation Methods"""
def calculateG(sequence,temperature):
    temperature += 273.15
    totalG = 0

    GInitiation = 0.2 - (temperature*-5.7/1000) # in kcal/mol

    GTermination = 0
    if sequence[len(sequence)-1] == "A" or sequence[len(sequence)-1] == "T":
        GTermination = 2.2 - (temperature*6.9/1000) # in kcal/mol

    totalG = totalG + GInitiation + GTermination

    for i in range(len(sequence)-1):
        bases = sequence[i:i+2]
        enthalpy = findH(bases) # in kcal/mol
        entropy = findS(bases) # in cal/(mol K)
        G = enthalpy - (temperature*entropy/1000) # in kcal/mol
        totalG += G

    return totalG

def calculateH(sequence):
	totalH = 0
	for i in range(len(sequence)-1):
		bases = sequence[i:i+2]
		enthalpy = findH(bases)
		totalH += enthalpy
	return totalH

def calculateS(sequence):
	totalS = 0
	for i in range(len(sequence)-1):
		bases = sequence[i:i+2]
		entropy = findS(bases)
		totalS += entropy
	return totalS

def findH(bases):
	if bases == "AA":
		return -7.6
	elif bases == "TT":
		return -7.6 # error
	elif bases == "AT":
		return -7.2
	elif bases == "TA":
		return -7.2
	elif bases == "CA":
		return -8.5
	elif bases == "AC":
		return -8.5 # error
	elif bases == "GT":
		return -8.4
	elif bases == "TG":
		return -8.4 # error
	elif bases == "CT":
		return -7.8
	elif bases == "TC":
		return -7.8 # error
	elif bases == "GA":
		return -8.2
	elif bases == "AG":
		return -8.2 # error
	elif bases == "CG":
		return -10.6
	elif bases == "GC":
		return -9.8
	elif bases == "GG":
		return -8.0
	elif bases == "CC":
		return -8.0 # error
	else:
		print("error")
		return 10000000

def findS(bases):
	if bases == "AA":
		return -21.3
	elif bases == "TT":
		return -21.3 # error
	elif bases == "AT":
		return -20.4
	elif bases == "TA":
		return -21.3
	elif bases == "CA":
		return -22.7
	elif bases == "AC":
		return -22.7 # error
	elif bases == "GT":
		return -22.4
	elif bases == "TG":
		return -22.4 # error
	elif bases == "CT":
		return -21.0
	elif bases == "TC":
		return -21.0 # error
	elif bases == "GA":
		return -22.2
	elif bases == "AG":
		return -22.2 # error
	elif bases == "CG":
		return -27.2
	elif bases == "GC":
		return -24.4
	elif bases == "GG":
		return -19.9
	elif bases == "CC":
		return -19.9 # error
	else:
		print("error")
		return 10000000

"""Graphs with ΔG"""
def TvsG(GValues,temperature):
    plt.plot(temperature,GValues[0],label='2 SE + 24 arms',color='khaki',marker='o',markersize=1)
    plt.plot(temperature,GValues[1],label='4 SE + 12 arms',color='gold',marker='o',markersize=1)
    plt.plot(temperature,GValues[2],label='6 SE + 8 arms',color='goldenrod',marker='o',markersize=1)
    plt.plot(temperature,GValues[3],label='8 SE + 6 arms',color='darkgoldenrod',marker='o',markersize=1)
    plt.plot(temperature,GValues[4],label='12 SE + 4 arms',color='saddlebrown',marker='o',markersize=1)
    plt.plot(temperature,GValues[5],label='16 SE + 3 arms',color='brown',marker='o',markersize=1)

    plt.xlabel("Temperature (°C)")
    plt.ylabel("ΔG (kcal/mol)")
    plt.title("Estimated ΔG values for corresponding transition temperatures")
    plt.legend(title='# of sticky ends')

    plt.show()

def GvsT(GValues,temperature):
    plt.plot(GValues[0],temperature,label='2 SE + 24 arms',color='khaki',marker='o',markersize=1)
    plt.plot(GValues[1],temperature,label='4 SE + 12 arms',color='gold',marker='o',markersize=1)
    plt.plot(GValues[2],temperature,label='6 SE + 8 arms',color='goldenrod',marker='o',markersize=1)
    plt.plot(GValues[3],temperature,label='8 SE + 6 arms',color='darkgoldenrod',marker='o',markersize=1)
    plt.plot(GValues[4],temperature,label='12 SE + 4 arms',color='saddlebrown',marker='o',markersize=1)
    plt.plot(GValues[5],temperature,label='16 SE + 3 arms',color='brown',marker='o',markersize=1)

    plt.xlabel("ΔG (kcal/mol)")
    plt.ylabel("Temperature (°C)",)
    plt.title("Estimated transition temperatures for corresponding ΔG values")
    plt.legend(title='# of sticky ends')

    plt.show()

def TvsGwithStickyNucleotides(GValues,temperature):
    plt.plot(temperature,GValues[0],label='2 SE + 24 arms',color='khaki',marker='o',markersize=1)
    plt.plot(temperature,GValues[1],label='4 SE + 12 arms',color='gold',marker='o',markersize=1)
    plt.plot(temperature,GValues[2],label='6 SE + 8 arms',color='goldenrod',marker='o',markersize=1)
    plt.plot(temperature,GValues[3],label='8 SE + 6 arms',color='darkgoldenrod',marker='o',markersize=1)
    plt.plot(temperature,GValues[4],label='12 SE + 4 arms',color='saddlebrown',marker='o',markersize=1)
    plt.plot(temperature,GValues[5],label='16 SE + 3 arms',color='brown',marker='o',markersize=1)

    plt.xlabel("Temperature (°C)")
    plt.ylabel("ΔG (kcal/mol)")
    plt.title("Estimated ΔG values of sticky nucleotides for corresponding transition temperatures")
    plt.legend(title='# of sticky ends + # of arms')

    plt.show()


"""Graphs with ΔH"""
def StickyEndvsH(HValues,sequence):
	stickyEnds = []
	for i in sequence:
		stickyEnds.append(len(i))

	plt.plot(stickyEnds,HValues,label="NN model",color='blue',marker='o',markersize=5,linestyle="None")

	plt.xlabel("# of sticky ends")
	plt.ylabel("ΔH (kcal/mol)")
	plt.title("Estimated ΔH values for corresponding # of sticky ends")

	plt.show()

def HwithStickyNucleotides(HValues,sequence):
	stickyEnds = []
	for i in sequence:
		stickyEnds.append(len(i))

	plt.plot(stickyEnds,HValues,label="NN model",color='blue',marker='o',markersize=5,linestyle="None")

	plt.xlabel("# of sticky ends")
	plt.ylabel("ΔH (kcal/mol)")
	plt.title("Estimated ΔH values for corresponding # of sticky ends taking # of arms into account")

	plt.show()

"""Graphs with ΔS"""
def StickyEndvsS(SValues,sequence):
	stickyEnds = []
	for i in sequence:
		stickyEnds.append(len(i))

	plt.plot(stickyEnds,SValues,label="NN model",color='springgreen',marker='o',markersize=5,linestyle="None")

	plt.xlabel("# of sticky ends")
	plt.ylabel("ΔS (cal mol-1 K-1)")
	plt.title("Estimated ΔS values for # of sticky ends")

	plt.show()

def SwithStickyNucleotides(SValues,sequence):
    stickyEnds = []
    for i in sequence:
    	stickyEnds.append(len(i))

    plt.plot(stickyEnds,SValues,label="NN model",color='springgreen',marker='o',markersize=5,linestyle="None")

    plt.xlabel("# of sticky ends")
    plt.ylabel("ΔS (cal mol-1 K-1)")
    plt.title("Estimated ΔS values for corresponding # of sticky ends taking # of arms into account")

    plt.show()


main()
