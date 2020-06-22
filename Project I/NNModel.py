import math
import numpy as np
import matplotlib.pyplot as plt

def main():
    sequence = ['GCGC','GCTAGC','GCTCGAGC','GACTCGAGTC','GCTAGCGCTAGC']

    TakinoueStickyEnd = [4,6,8,10,12]
    TakinoueH = [-30.0,-42.2,-62.0,-75.6,-95.2]
    TakinoueT = [46.3,48.3,63.7,61.3,70.3]

    temperature = np.arange(20,80,0.1).tolist()

    GValues = []
    for s in sequence:
        GValue = []
        for t in temperature:
            GValue.append(calculateG(s,t))
        GValues.append(GValue)

    TakinoueG = [calculateG(sequence[0],TakinoueT[0]),calculateG(sequence[1],TakinoueT[1]),calculateG(sequence[2],TakinoueT[2]),calculateG(sequence[3],TakinoueT[3]),calculateG(sequence[4],TakinoueT[4])]

    HValues = [calculateH(sequence[0]),calculateH(sequence[1]),calculateH(sequence[2]),calculateH(sequence[3]),calculateH(sequence[4])]
    SValues = [calculateS(sequence[0]),calculateS(sequence[1]),calculateS(sequence[2]),calculateS(sequence[3]),calculateS(sequence[4])]


    """----------------------------ΔG----------------------------"""
    """Estimated ΔG values for corresponding transition temperatures"""
    TvsG(GValues,temperature)

    """Estimated ΔG values for corresponding transition temperatures with Takinoue's (2020) data"""
    #TvsGwithTakinoueG(GValues,temperature,TakinoueG,TakinoueT)

    """Estimated transition temperatures for corresponding ΔG values"""
    #GvsT(GValues,temperature)

    """Estimated transition temperatures for corresponding ΔG values with Takinoue's (2020) data"""
    #GvsTwithTakinoueG(GValues,temperature,TakinoueG,TakinoueT)


    """----------------------------ΔH----------------------------"""
    """Estimated ΔH values for corresponding # of sticky ends"""
    #StickyEndvsH(HValues,sequence)

    """Comparing the relationship between sticky and ΔH with Takinoue's (2020) data"""
    #StickyEndvsHwithTakinoueH(HValues,sequence,TakinoueH)

    """Estimated transition temperatures for corresponding ΔH"""
    """Comparing closeness of data with Takinoue (2020)"""
    #HvsTwithTakinoueH(HValues,TakinoueT,TakinoueH)

    """Alternative graph with -ΔH"""
    #NegativeHvsTwithTakinoueH(HValues,TakinoueT,TakinoueH)


    """----------------------------ΔS----------------------------"""
    """Estimated ΔS values for corresponding # of sticky ends"""
    #StickyEndvsS(SValues,sequence)

    """Estimated transition temperatures for corresponding ΔS"""
    #SvsTwithTakinoueT(SValues,TakinoueT)


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
    plt.plot(temperature,GValues[0],label='4',color='lemonchiffon',marker='o',markersize=1)
    plt.plot(temperature,GValues[1],label='6',color='khaki',marker='o',markersize=1)
    plt.plot(temperature,GValues[2],label='8',color='gold',marker='o',markersize=1)
    plt.plot(temperature,GValues[3],label='10',color='goldenrod',marker='o',markersize=1)
    plt.plot(temperature,GValues[4],label='12',color='darkgoldenrod',marker='o',markersize=1)

    plt.xlabel("Temperature (°C)")
    plt.ylabel("ΔG (kcal/mol)")
    plt.title("Estimated ΔG values for corresponding transition temperatures")
    plt.legend(title='# of sticky ends')

    plt.show()

def TvsGwithTakinoueG(GValues,temperature,TakinoueG,TakinoueT):
    plt.plot(temperature,GValues[0],label='4',color='lemonchiffon',marker='o',markersize=1)
    plt.plot(temperature,GValues[1],label='6',color='khaki',marker='o',markersize=1)
    plt.plot(temperature,GValues[2],label='8',color='gold',marker='o',markersize=1)
    plt.plot(temperature,GValues[3],label='10',color='goldenrod',marker='o',markersize=1)
    plt.plot(temperature,GValues[4],label='12',color='darkgoldenrod',marker='o',markersize=1)
    plt.plot(TakinoueT,TakinoueG,color='orange',marker='o',markersize=5,linestyle='None')

    plt.xlabel("Temperature (°C)")
    plt.ylabel("ΔG (kcal/mol)")
    plt.title("Estimated ΔG values for corresponding transition temperatures with Takinoue's data")
    plt.legend(title='# of sticky ends')

    plt.show()

def GvsT(GValues,temperature):
    plt.plot(GValues[0],temperature,label='4',color='lemonchiffon',marker='o',markersize=1)
    plt.plot(GValues[1],temperature,label='6',color='khaki',marker='o',markersize=1)
    plt.plot(GValues[2],temperature,label='8',color='gold',marker='o',markersize=1)
    plt.plot(GValues[3],temperature,label='10',color='goldenrod',marker='o',markersize=1)
    plt.plot(GValues[4],temperature,label='12',color='darkgoldenrod',marker='o',markersize=1)

    plt.xlabel("ΔG (kcal/mol)")
    plt.ylabel("Temperature (°C)",)
    plt.title("Estimated transition temperatures for corresponding ΔG values")
    plt.legend(title='# of sticky ends')

    plt.show()

def GvsTwithTakinoueG(GValues,temperature,TakinoueG,TakinoueT):
    plt.plot(GValues[0],temperature,label='4',color='lemonchiffon',marker='o',markersize=1)
    plt.plot(GValues[1],temperature,label='6',color='khaki',marker='o',markersize=1)
    plt.plot(GValues[2],temperature,label='8',color='gold',marker='o',markersize=1)
    plt.plot(GValues[3],temperature,label='10',color='goldenrod',marker='o',markersize=1)
    plt.plot(GValues[4],temperature,label='12',color='darkgoldenrod',marker='o',markersize=1)
    plt.plot(TakinoueG,TakinoueT,color='orange',marker='o',markersize=5,linestyle='None')

    plt.xlabel("ΔG (kcal/mol)")
    plt.ylabel("Temperature (°C)",)
    plt.title("Estimated transition temperatures for corresponding ΔG values with Takinoue's data")
    plt.legend(title='# of sticky ends')

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

def StickyEndvsHwithTakinoueH(HValues,sequence,TakinoueH):
	stickyEnds = []
	for i in sequence:
		stickyEnds.append(len(i))

	plt.plot(stickyEnds,HValues,label="NN model",color='blue',marker='o',markersize=5,linestyle="None")
	plt.plot(stickyEnds,TakinoueH,label="Takinoue (2020)",color="cyan",marker='o',markersize=5,linestyle="None")

	plt.xlabel("# of sticky ends")
	plt.ylabel("ΔH (kcal/mol)")
	plt.title("Estimated ΔH values for corresponding # of sticky ends with Takinoue's (2020) data")
	plt.legend()

	plt.show()

def HvsTwithTakinoueH(HValues,TakinoueT,TakinoueH):
	plt.plot(HValues,TakinoueT,label="NN model",color='blue',marker='o',markersize=5,linestyle="None")
	plt.plot(TakinoueH,TakinoueT,label="Takinoue (2020)",color="cyan",marker='o',markersize=5,linestyle="None")

	plt.xlabel("ΔH (kcal/mol)")
	plt.ylabel("Temperature (°C)")
	plt.title("Estimated transition temperatures for corresponding ΔH")
	plt.legend()

	plt.show()

def NegativeHvsTwithTakinoueH(HValues,TakinoueT,TakinoueH):
    for i in HValues:
        i *= -1

    for i in TakinoueH:
        i *= -1

    plt.plot(HValues,TakinoueT,label="NN model",color='blue',marker='o',markersize=5,linestyle="None")
    plt.plot(TakinoueH,TakinoueT,label="Takinoue (2020)",color="cyan",marker='o',markersize=5,linestyle="None")

    plt.xlabel("-ΔH (kcal/mol)")
    plt.ylabel("Temperature (°C)")
    plt.title("Estimated transition temperatures for corresponding ΔH")
    plt.legend()

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

def SvsTwithTakinoueT(SValues,TakinoueT):
	plt.plot(SValues,TakinoueT,color='springgreen',marker='o',markersize=5,linestyle="None")

	plt.xlabel("ΔS (cal mol-1 K-1)")
	plt.ylabel("Temperature (°C)")
	plt.title("Estimated transition temperatures for corresponding ΔS")

	plt.show()


main()
