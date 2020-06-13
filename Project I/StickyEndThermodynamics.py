"""calculates ∆G from transition temperatures using nearest-neighbor thermodynamic parameters for DNA Watson-Crick pairs in 1 M NaCl"""

def main():
	sequence = str(input("Sequence (5' to 3'): "))
	temperature = 273.15 + float(input("Temperature (C): "))

	totalGibbs = 0

	GibbsInitiation = 0.2 - (temperature*-5.7/1000) # in kcal/mol

	GibbsTermination = 0
	if sequence[len(sequence)-1] == "A" or sequence[len(sequence)-1] == "T":
		GibbsTermination = 2.2 - (temperature*6.9/1000) # in kcal/mol

	totalGibbs = totalGibbs + GibbsInitiation + GibbsTermination

    # checkBase = ""

	for i in range(len(sequence)-1):
		bases = sequence[i:i+2]
		enthalpy = findH(bases) # in kcal/mol
		entropy = findS(bases) # in cal/(mol K)
		Gibbs = enthalpy - (temperature*entropy/1000) # in kcal/mol
		totalGibbs += Gibbs

	roundedGibbs = round(totalGibbs,3)

	print(roundedGibbs, "kcal/mol")
	# print("Note: did not take into account symmetry correction some literature values were not included and were thus estimated.")

# taken from THE THERMODYNAMICS OF DNA STRUCTURAL MOTIFS (2004) by SantaLucia & Hicks
# some values marked "not sure" were estimated because they were not included in the literature
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

# taken from THE THERMODYNAMICS OF DNA STRUCTURAL MOTIFS (2004) by SantaLucia & Hicks
# some values marked "not sure" were estimated because they were not included in the literature
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

main()
