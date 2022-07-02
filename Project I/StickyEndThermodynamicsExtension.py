import math
import numpy as np
import matplotlib.pyplot as plt

def main():
	# sequence = str(input("Sequence (5' to 3'): "))
	sequence = ["GCTAGCGCTAGC","GACTCGAGTC","GCTCGAGC","GCTAGC","GCGC"]

	temperature = np.arange(20, 80, 0.1).tolist()

	dataTakinoueStickyEndvsH = [[4,6,8,10,12], [-30.0,-42.2,-62.0,-75.6,-95.2]]
	dataTakinoueStickyEndvsT = [[4,6,8,10,12], [46.3,48.3,63.7,61.3,70.3]]

	dataMelissaArms = [[4,6,8,12,24], [29,43,46,52,55]] # for CGATCG
	dataTakinoueArms = [[3,4,6], [65, 69, 75]] # for GCTCGAGC

	"""----------------------------ΔG----------------------------"""
	GValues = [findGValues(sequence[0], temperature), findGValues(sequence[1], temperature), findGValues(sequence[2], temperature), findGValues(sequence[3], temperature), findGValues(sequence[4], temperature)]
	#GvsT(GValues, temperature)
	#TvsG(GValues, temperature)
	#GvsTWithData(GValues, temperature, dataMelissaArms, dataTakinoueArms)

	"""ΔG for Takinoue (2019)"""
	TakinoueValues = [-77.75715000000001, -62.57411999999999, -51.29203, -36.8523, -26.21729]
	#xyPlot(TakinoueValues, [70.3,61.3,63.7,48.3,46.3])
	"""----------------------------ΔG----------------------------"""


	"""----------------------------ΔH----------------------------"""
	HValues = [findHValues(sequence[0]), findHValues(sequence[1]), findHValues(sequence[2]), findHValues(sequence[3]), findHValues(sequence[4])]
	"""Relationship between sticky end and ΔH"""
	#StickyEndvsH(HValues, sequence)

	"""Comparing the relationship between sticky and ΔH with Takinoue's (2019) data"""
	#StickyEndvsHWithData(HValues, sequence, dataTakinoueStickyEndvsH)

	"""Figure 1: Estimated transition temperatures for corresponding ΔH + Comparing closeness of data with Takinoue (2019)"""
	dataProgramHvsT = [[HValues],[46.3,48.3,63.7,61.3,70.3]]
	dataTakinoueHvsT = [[-30.0,-42.2,-62.0,-75.6,-95.2],[46.3,48.3,63.7,61.3,70.3]]
	#HvsTwithData(dataProgramHvsT,dataTakinoueHvsT)

	"""Alternative Figure 1 with -ΔH"""
	dataProgramNegativeHvsT = [[30.2,42.8,62.2,75.5,96.2],[46.3,48.3,63.7,61.3,70.3]]
	dataTakinoueNegativeHvsT = [[30.0,42.2,62.0,75.6,95.2],[46.3,48.3,63.7,61.3,70.3]]
	#NegativeHvsTwithData(dataProgramNegativeHvsT,dataTakinoueNegativeHvsT)
	"""----------------------------ΔH----------------------------"""


	"""----------------------------ΔS----------------------------"""
	SValues = [findSValues(sequence[0]), findSValues(sequence[1]), findSValues(sequence[2]), findSValues(sequence[3]), findSValues(sequence[4])]
	"""Relationship between sticky end and ΔS"""
	#StickyEndvsS(SValues, sequence)

	"""Figure 2: Estimated transition temperatures for corresponding ΔS"""
	dataProgramSvsT = [SValues,[46.3,48.3,63.7,61.3,70.3]]
	#SvsTwithData(dataProgramSvsT)
	"""----------------------------ΔS----------------------------"""


def findGValues(sequence, temperature):
	GValues = []
	for temp in temperature:
		temp += 273.15

		totalGibbs = 0

		GibbsInitiation = 0.2 - (temp*-5.7/1000) # in kcal/mol

		GibbsTermination = 0
		if sequence[len(sequence)-1] == "A" or sequence[len(sequence)-1] == "T":
			GibbsTermination = 2.2 - (temp*6.9/1000) # in kcal/mol

		totalGibbs = totalGibbs + GibbsInitiation + GibbsTermination

		# checkBase = ""

		for i in range(len(sequence)-1):
			bases = sequence[i:i+2]
			enthalpy = findH(bases) # in kcal/mol
			entropy = findS(bases) # in cal/(mol K)
			Gibbs = enthalpy - (temp*entropy/1000) # in kcal/mol
			totalGibbs += Gibbs

		# roundedGibbs = round(totalGibbs, 3-int(math.floor(math.log10(abs(totalGibbs))))-1)

		GValues.append(totalGibbs)

	return GValues

def findHValues(sequence):
	totalH = 0
	for i in range(len(sequence)-1):
		bases = sequence[i:i+2]
		enthalpy = findH(bases)
		totalH += enthalpy
	return totalH

def findSValues(sequence):
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
		return -7.6
        # checkBase += "*"
	elif bases == "AT":
		return -7.2
	elif bases == "TA":
		return -7.2
	elif bases == "CA":
		return -8.5
	elif bases == "AC":
		return -8.5
        # checkBase += "*"
	elif bases == "GT":
		return -8.4
	elif bases == "TG":
		return -8.4
        # checkBase += "*"
	elif bases == "CT":
		return -7.8
	elif bases == "TC":
		return -7.8
        # checkBase += "*"
	elif bases == "GA":
		return -8.2
	elif bases == "AG":
		return -8.2
        # checkBase += "*"
	elif bases == "CG":
		return -10.6
	elif bases == "GC":
		return -9.8
	elif bases == "GG":
		return -8.0
	elif bases == "CC":
		return -8.0
        # checkBase += "*"
	else:
		print("error")
		return 10000000

def findS(bases):
	if bases == "AA":
		return -21.3
	elif bases == "TT":
		return -21.3
        # checkBase += "*"
	elif bases == "AT":
		return -20.4
	elif bases == "TA":
		return -21.3
	elif bases == "CA":
		return -22.7
	elif bases == "AC":
		return -22.7
        # checkBase += "*"
	elif bases == "GT":
		return -22.4
	elif bases == "TG":
		return -22.4
        # checkBase += "*"
	elif bases == "CT":
		return -21.0
	elif bases == "TC":
		return -21.0
        # checkBase += "*"
	elif bases == "GA":
		return -22.2
	elif bases == "AG":
		return -22.2
        # checkBase += "*"
	elif bases == "CG":
		return -27.2
	elif bases == "GC":
		return -24.4
	elif bases == "GG":
		return -19.9
	elif bases == "CC":
		return -19.9
        # checkBase += "*"
	else:
		print("error")
		return 10000000

"""Plots with G"""
def GvsT(GValues, temperature):
	plt.plot(GValues[0], temperature, label='12', color='blue', marker='o', markersize=1)
	plt.plot(GValues[1], temperature, label='10', color='green', marker='o', markersize=1)
	plt.plot(GValues[2], temperature, label='8', color='red', marker='o', markersize=1)
	plt.plot(GValues[3], temperature, label='6', color='cyan', marker='o', markersize=1)
	plt.plot(GValues[4], temperature, label='4', color='magenta', marker='o', markersize=1)

	plt.xlabel("ΔG (kcal/mol)")
	plt.ylabel("Temperature (°C)")
	plt.title("Estimated transition temperatures for ΔG values")
	plt.legend(title='# of sticky ends')

	plt.show()

def TvsG(GValues, temperature):
	plt.plot(temperature, GValues[0], label='12', color='blue', marker='o', markersize=1)
	plt.plot(temperature, GValues[1], label='10', color='green', marker='o', markersize=1)
	plt.plot(temperature, GValues[2], label='8', color='red', marker='o', markersize=1)
	plt.plot(temperature, GValues[3], label='6', color='cyan', marker='o', markersize=1)
	plt.plot(temperature, GValues[4], label='4', color='magenta', marker='o', markersize=1)

	plt.xlabel("Temperature (°C)")
	plt.ylabel("ΔG (kcal/mol)")
	plt.title("Estimated ΔG values for transition temperatures")
	plt.legend(title='# of sticky ends')

	plt.show()

def GvsTWithData(GValues, temperature, data1, data2):
	plt.plot(GValues[0], temperature, label='12', color='blue', marker='o', markersize=1)
	plt.plot(GValues[1], temperature, label='10', color='green', marker='o', markersize=1)
	plt.plot(GValues[2], temperature, label='8', color='red', marker='o', markersize=1)
	plt.plot(GValues[3], temperature, label='6', color='cyan', marker='o', markersize=1)
	plt.plot(GValues[4], temperature, label='4', color='magenta', marker='o', markersize=1)

	plt.xlabel("ΔG (kcal/mol)")
	plt.ylabel("Temperature (°C)")
	plt.title("Estimated transition temperatures for ΔG values")
	plt.legend(title='# of sticky ends')


	plt.show()

def xyPlot(GValues, temperature):
	plt.plot(GValues, temperature, color='blue', marker='o', markersize=5, linestyle="None")
	plt.xlabel("ΔG (kcal/mol)")
	plt.ylabel("Temperature (°C)")
	plt.title("Estimated transition temperatures for corresponding ΔG")

	z = np.polyfit(GValues, temperature, 1)
	p = np.poly1d(z)
	plt.plot(GValues,p(GValues),"b--")
	print("Figure 3: y=%.6fx+(%.6f)"%(z[0],z[1]))

	plt.show()

"""Plots with H"""
def StickyEndvsH(HValues, sequence):
	stickyEnds = []
	for i in sequence:
		stickyEnds.append(len(i))

	plt.plot(stickyEnds, HValues, label="program", color='blue', marker='o', markersize=5, linestyle="None")

	plt.xlabel("# of sticky ends")
	plt.ylabel("ΔH (kcal/mol)")
	plt.title("Estimated ΔH values for # of sticky ends")

	z = np.polyfit(stickyEnds, HValues, 1)
	p = np.poly1d(z)
	plt.plot(stickyEnds,p(stickyEnds),"r--")
	print("y=%.6fx+(%.6f)"%(z[0],z[1]))

	plt.show()

def StickyEndvsHWithData(HValues, sequence, data):
	stickyEnds = []
	for i in sequence:
		stickyEnds.append(len(i))

	plt.plot(stickyEnds, HValues, label="program", color='blue', marker='o', markersize=5, linestyle="None")
	plt.plot(data[0], data[1], label="Takinoue (2019)", color="cyan", marker='o', markersize=5, linestyle="None")

	plt.xlabel("# of sticky ends")
	plt.ylabel("ΔH (kcal/mol)")
	plt.title("Estimated ΔH values for # of sticky ends")
	plt.legend()

	z1 = np.polyfit(stickyEnds, HValues, 1)
	p1 = np.poly1d(z1)
	plt.plot(stickyEnds,p1(stickyEnds),"b--")
	print("Program: y=%.6fx+(%.6f)"%(z1[0],z1[1]))

	z2 = np.polyfit(data[0], data[1], 1)
	p2 = np.poly1d(z2)
	plt.plot(data[0],p2(data[0]),"c--")
	print("Takinoue (2019): y=%.6fx+(%.6f)"%(z2[0],z2[1]))

	plt.show()

def HvsTwithData(data1,data2):
	plt.plot(data1[0], data1[1], label="program", color='blue', marker='o', markersize=5, linestyle="None")
	plt.plot(data2[0], data2[1], label="Takinoue (2019)", color="cyan", marker='o', markersize=5, linestyle="None")

	plt.xlabel("ΔH (kcal/mol)")
	plt.ylabel("Temperature (°C)")
	plt.title("Estimated transition temperatures for corresponding ΔH")
	plt.legend()

	z1 = np.polyfit(data1[0], data1[1], 1)
	p1 = np.poly1d(z1)
	plt.plot(data1[0],p1(data1[0]),"b--")
	print("Figure 1: Program  y=%.6fx+(%.6f)"%(z1[0],z1[1]))

	z2 = np.polyfit(data2[0], data2[1], 1)
	p2 = np.poly1d(z2)
	plt.plot(data2[0],p2(data2[0]),"c--")
	print("Figure 1: Takinoue (2019)  y=%.6fx+(%.6f)"%(z2[0],z2[1]))

	plt.show()

def NegativeHvsTwithData(data1,data2):
	plt.plot(data1[0], data1[1], label="program", color='blue', marker='o', markersize=5, linestyle="None")
	plt.plot(data2[0], data2[1], label="Takinoue (2019)", color="cyan", marker='o', markersize=5, linestyle="None")

	plt.xlabel("-ΔH (kcal/mol)")
	plt.ylabel("Temperature (°C)")
	plt.title("Estimated transition temperatures for corresponding ΔH")
	plt.legend()

	z1 = np.polyfit(data1[0], data1[1], 1)
	p1 = np.poly1d(z1)
	plt.plot(data1[0],p1(data1[0]),"b--")
	print("Program: y=%.6fx+(%.6f)"%(z1[0],z1[1]))

	z2 = np.polyfit(data2[0], data2[1], 1)
	p2 = np.poly1d(z2)
	plt.plot(data2[0],p2(data2[0]),"c--")
	print("Takinoue (2019): y=%.6fx+(%.6f)"%(z2[0],z2[1]))

	plt.show()

"""Plots with S"""
def StickyEndvsS(SValues, sequence):
	stickyEnds = []
	for i in sequence:
		stickyEnds.append(len(i))

	plt.plot(stickyEnds, SValues, label="program", color='blue', marker='o', markersize=5, linestyle="None")

	plt.xlabel("# of sticky ends")
	plt.ylabel("ΔS (cal mol-1 K-1)")
	plt.title("Estimated ΔS values for # of sticky ends")

	z = np.polyfit(stickyEnds, SValues, 1)
	p = np.poly1d(z)
	plt.plot(stickyEnds,p(stickyEnds),"r--")
	print("y=%.6fx+(%.6f)"%(z[0],z[1]))

	plt.show()

def SvsTwithData(data):
	plt.plot(data[0], data[1], color='blue', marker='o', markersize=5, linestyle="None")

	plt.xlabel("ΔS (cal mol-1 K-1)")
	plt.ylabel("Temperature (°C)")
	plt.title("Estimated transition temperatures for corresponding ΔS")

	z = np.polyfit(data[0], data[1], 1)
	p = np.poly1d(z)
	plt.plot(data[0],p(data[0]),"b--")
	print("Figure 2: y=%.6fx+(%.6f)"%(z[0],z[1]))

	plt.show()


main()


def main2():
	Msequence1 = "CGATCG"
	Msequence2 = "GCGC" #4 arms, 26 °C
	Tsequence1 = "GCTCGAGC"
	Tsequence2 = ["GCTAGCGCTAGC","GACTCGAGTC","GCTCGAGC","GCTAGC","GCGC"]

	# 6-nucleotide sequence
	MelissaArms = [4,6,8,12,24]
	MelissaT = [29,43,46,52,55]

	TakinoueArms = [3,4,6]
	TakinoueT = [65,69,75]

	MelissaSNt = [24,36,48,72,144]
	TakinoueSNt = [24,32,48]
	TakinoueSNt2 = [12,18,24,30,36]
	TakinoueT2 = [46.3,48.3,63.7,61.3,70.3]

	#ArmsvsT([MelissaArms, MelissaT],[TakinoueArms, TakinoueT])
	#SNtvsT([MelissaSNt, MelissaT],[TakinoueSNt, TakinoueT],[TakinoueSNt2, TakinoueT2])

	"""ΔH"""
	Msequence1H = -44.4
	Msequence2H = -30.2
	Tsequence1H = -62.2
	Tsequence2H = [-30.2, -42.8, -62.2, -75.5, -96.2]

	Msequence1HSNt = [-177.6, -266.4, -355.2, -532.8, -1065.6]
	Tsequence1HSNt = [-186.6, -248.8, -373.2]
	Tsequence2HSNt = [-90.6, -128.4, -186.6, -226.5, -288.6]

	#HSNtvsT([Msequence1HSNt,MelissaT],[Tsequence1HSNt,TakinoueT],[Tsequence2HSNt,TakinoueT2])

	Msequence1NegativeHSNt = [177.6, 266.4, 355.2, 532.8, 1065.6]
	Tsequence1NegativeHSNt = [186.6, 248.8, 373.2]
	Tsequence2NegativeHSNt = [90.6, 128.4, 186.6, 226.5, 288.6]

	#NegativeHSNtvsT([Msequence1NegativeHSNt,MelissaT],[Tsequence1NegativeHSNt,TakinoueT],[Tsequence2NegativeHSNt,TakinoueT2])

	"""ΔS"""
	Msequence1S = -118.0
	Msequence2S = -76.0
	Tsequence1S = -162.4
	Tsequence2S = [-76.0, -113.3, -162.4, -201.9, -253.8]

	Msequence1SSNt = [-472.0, -708.0, -944.0, -1416.0, -2832.0]
	Tsequence1SSNt = [-487.2, -649.6, -974.4]
	Tsequence2SSNt = [-228.0, -339.9, -487.2, -605.7, -761.4]

	SSNtvsT([Msequence1SSNt,MelissaT],[Tsequence1SSNt,TakinoueT],[Tsequence2SSNt,TakinoueT2])


	"""ΔG"""
	##Msequence1G = Msequence1H -
	#Msequence2G =


def ArmsvsT(data1,data2):
	plt.plot(data1[0], data1[1], label="Melissa's data (6 nt)", color='blue', marker='o', markersize=5, linestyle="None")
	plt.plot(data2[0], data2[1], label="Takinoue's data (8 nt)", color="cyan", marker='o', markersize=5, linestyle="None")

	plt.xlabel("# of arms")
	plt.ylabel("Temperature (°C)")
	plt.title("Estimated transition temperatures for corresponding # of arms")
	plt.legend()

	z1 = np.polyfit(data1[0], data1[1], 1)
	p1 = np.poly1d(z1)
	plt.plot(data1[0],p1(data1[0]),"b--")
	print("Melissa's data: y=%.6fx+(%.6f)"%(z1[0],z1[1]))

	z2 = np.polyfit(data2[0], data2[1], 1)
	p2 = np.poly1d(z2)
	plt.plot(data2[0],p2(data2[0]),"c--")
	print("Takinoue's data: y=%.6fx+(%.6f)"%(z2[0],z2[1]))

	plt.show()

def SNtvsT(data1,data2,data3):
	plt.plot(data1[0], data1[1], label="Melissa's data (4,6,8,12,24 arms: 6 nt)", color='blue', marker='o', markersize=5, linestyle="None")
	plt.plot(16, 26, label="Melissa's data (4 arms: 4 nt)", color="red", marker='o', markersize=5, linestyle="None")
	plt.plot(data2[0], data2[1], label="Takinoue's data (3,4,6 arms: 8 nt)", color="cyan", marker='o', markersize=5, linestyle="None")
	plt.plot(data3[0], data3[1], label="Takinoue's data (3 arms: 4,6,8,10,12 nt)", color="green", marker='o', markersize=5, linestyle="None")

	plt.xlabel("# of 'sticky nucleotides'")
	plt.ylabel("Temperature (°C)")
	plt.title("Estimated transition temperatures for corresponding # of 'sticky nucleotides'")
	plt.legend()

	z1 = np.polyfit(data1[0], data1[1], 1)
	p1 = np.poly1d(z1)
	plt.plot(data1[0],p1(data1[0]),"b--")
	print("Melissa's data: y=%.6fx+(%.6f)"%(z1[0],z1[1]))

	z2 = np.polyfit(data2[0], data2[1], 1)
	p2 = np.poly1d(z2)
	plt.plot(data2[0],p2(data2[0]),"c--")
	print("Takinoue's data: y=%.6fx+(%.6f)"%(z2[0],z2[1]))

	z3 = np.polyfit(data3[0], data3[1], 1)
	p3 = np.poly1d(z3)
	plt.plot(data3[0],p3(data3[0]),"g--")
	print("Takinoue's data: y=%.6fx+(%.6f)"%(z3[0],z3[1]))

	plt.show()

def HSNtvsT(data1,data2,data3):
	plt.plot(data1[0], data1[1], label="Melissa's data (4,6,8,12,24 arms: 6 nt)", color='blue', marker='o', markersize=5, linestyle="None")
	plt.plot(-120.8, 26, label="Melissa's data (4 arms: 4 nt)", color="red", marker='o', markersize=5, linestyle="None")
	plt.plot(data2[0], data2[1], label="Takinoue's data (3,4,6 arms: 8 nt)", color="cyan", marker='o', markersize=5, linestyle="None")
	plt.plot(data3[0], data3[1], label="Takinoue's data (3 arms: 4,6,8,10,12 nt)", color="green", marker='o', markersize=5, linestyle="None")

	plt.xlabel("ΔH (kcal/mol)")
	plt.ylabel("Temperature (°C)")
	plt.title("Estimated transition temperatures for corresponding ΔH (kcal/mol) multiplied by # of arms")
	plt.legend()

	z1 = np.polyfit(data1[0], data1[1], 1)
	p1 = np.poly1d(z1)
	plt.plot(data1[0],p1(data1[0]),"b--")
	print("Melissa's data: y=%.6fx+(%.6f)"%(z1[0],z1[1]))

	z2 = np.polyfit(data2[0], data2[1], 1)
	p2 = np.poly1d(z2)
	plt.plot(data2[0],p2(data2[0]),"c--")
	print("Takinoue's data: y=%.6fx+(%.6f)"%(z2[0],z2[1]))

	z3 = np.polyfit(data3[0], data3[1], 1)
	p3 = np.poly1d(z3)
	plt.plot(data3[0],p3(data3[0]),"g--")
	print("Takinoue's data: y=%.6fx+(%.6f)"%(z3[0],z3[1]))

	plt.show()

def NegativeHSNtvsT(data1,data2,data3):
	plt.plot(data1[0], data1[1], label="Melissa's data (4,6,8,12,24 arms: 6 nt)", color='blue', marker='o', markersize=5, linestyle="None")
	plt.plot(120.8, 26, label="Melissa's data (4 arms: 4 nt)", color="red", marker='o', markersize=5, linestyle="None")
	plt.plot(data2[0], data2[1], label="Takinoue's data (3,4,6 arms: 8 nt)", color="cyan", marker='o', markersize=5, linestyle="None")
	plt.plot(data3[0], data3[1], label="Takinoue's data (3 arms: 4,6,8,10,12 nt)", color="green", marker='o', markersize=5, linestyle="None")

	plt.xlabel("-ΔH (kcal/mol)")
	plt.ylabel("Temperature (°C)")
	plt.title("Estimated transition temperatures for corresponding ΔH (kcal/mol) multiplied by # of arms")
	plt.legend()

	z1 = np.polyfit(data1[0], data1[1], 1)
	p1 = np.poly1d(z1)
	plt.plot(data1[0],p1(data1[0]),"b--")
	print("Melissa's data: y=%.6fx+(%.6f)"%(z1[0],z1[1]))

	z2 = np.polyfit(data2[0], data2[1], 1)
	p2 = np.poly1d(z2)
	plt.plot(data2[0],p2(data2[0]),"c--")
	print("Takinoue's data: y=%.6fx+(%.6f)"%(z2[0],z2[1]))

	z3 = np.polyfit(data3[0], data3[1], 1)
	p3 = np.poly1d(z3)
	plt.plot(data3[0],p3(data3[0]),"g--")
	print("Takinoue's data: y=%.6fx+(%.6f)"%(z3[0],z3[1]))

	plt.show()

def SSNtvsT(data1,data2,data3):
	plt.plot(data1[0], data1[1], label="Melissa's data (4,6,8,12,24 arms: 6 nt)", color='blue', marker='o', markersize=5, linestyle="None")
	plt.plot(-304.0, 26, label="Melissa's data (4 arms: 4 nt)", color="red", marker='o', markersize=5, linestyle="None")
	plt.plot(data2[0], data2[1], label="Takinoue's data (3,4,6 arms: 8 nt)", color="cyan", marker='o', markersize=5, linestyle="None")
	plt.plot(data3[0], data3[1], label="Takinoue's data (3 arms: 4,6,8,10,12 nt)", color="green", marker='o', markersize=5, linestyle="None")

	plt.xlabel("ΔS (cal mol-1 K-1)")
	plt.ylabel("Temperature (°C)")
	plt.title("Estimated transition temperatures for corresponding ΔS (cal mol-1 K-1) multiplied by # of arms")
	plt.legend()

	z1 = np.polyfit(data1[0], data1[1], 1)
	p1 = np.poly1d(z1)
	plt.plot(data1[0],p1(data1[0]),"b--")
	print("Melissa's data: y=%.6fx+(%.6f)"%(z1[0],z1[1]))

	z2 = np.polyfit(data2[0], data2[1], 1)
	p2 = np.poly1d(z2)
	plt.plot(data2[0],p2(data2[0]),"c--")
	print("Takinoue's data: y=%.6fx+(%.6f)"%(z2[0],z2[1]))

	z3 = np.polyfit(data3[0], data3[1], 1)
	p3 = np.poly1d(z3)
	plt.plot(data3[0],p3(data3[0]),"g--")
	print("Takinoue's data: y=%.6fx+(%.6f)"%(z3[0],z3[1]))

	plt.show()

main2()
