import math
import numpy as np
import matplotlib.pyplot as plt

def main():
    MelissaSequence6 = "CGATCG"
    MelissaSequence4 = "GCGC"
    TakinoueSequence8 = "GCTCGAGC"
    TakinoueSequences = ['GCGC','GCTAGC','GCTCGAGC','GACTCGAGTC','GCTAGCGCTAGC']

    # 6-nucleotide sequence
    MelissaArms6 = [4,6,8,12,24]
    MelissaT6 = [29,43,46,52,55]

    # 4-nucleotide sequence
    MelissaArms4 = 4
    MelissaT4 = 26

    # 8-nucleotide sequence
    TakinoueArms8 = [3,4,6]
    TakinoueT8 = [65,69,75]

    # 4,6,8,10,12-nucleotide sequence
    TakinoueT = [46.3,48.3,63.7,61.3,70.3]

    # Sticky nucleotides = # of arms * # of sticky ends
    MelissaStickyNucleotides6 = [24,36,48,72,144]
    MelissaStickyNucleotides4 = 16
    TakinoueStickyNucleotides8 = [24,32,48]
    TakinoueStickyNucleotides = [12,18,24,30,36]


    """Estimated transition temperatures for corresponding # of 'sticky nucleotides'"""
    #StickyNucleotidesvsT([MelissaStickyNucleotides6,MelissaT6],[MelissaStickyNucleotides4,MelissaT4],[TakinoueStickyNucleotides8,TakinoueT8],[TakinoueStickyNucleotides,TakinoueT])


    """----------------------------ΔH----------------------------"""
    MelissaSequence6H = -44.4
    MelissaSequence4H = -30.2
    TakinoueSequence8H = -62.2
    TakinoueSequencesH = [-30.2,-42.8,-62.2,-75.5,-96.2]

    MelissaHStickyNucleotides6 = [-177.6,-266.4,-355.2,-532.8,-1065.6]
    MelissaHStickyNucleotides4 = -120.8
    TakinoueHStickyNucleotides8 = [-186.6,-248.8,-373.2]
    TakinoueHStickyNucleotides = [-90.6,-128.4,-186.6,-226.5,-288.6]

    """Estimated transition temperatures for corresponding ΔH multiplied by # of arms"""
    #HStickyNucleotidesvsT([MelissaHStickyNucleotides6,MelissaT6],[MelissaHStickyNucleotides4,MelissaT4],[TakinoueHStickyNucleotides8,TakinoueT8],[TakinoueHStickyNucleotides,TakinoueT])

    MelissaNegativeHStickyNucleotides6 = [177.6,266.4,355.2,532.8,1065.6]
    MelissaNegativeHStickyNucleotides4 = 120.8
    TakinoueNegativeHStickyNucleotides8 = [186.6,248.8,373.2]
    TakinoueNegativeHStickyNucleotides = [90.6,128.4,186.6,226.5,288.6]

    """Alternative graph with -ΔH"""
    #NegativeHStickyNucleotidesvsT([MelissaNegativeHStickyNucleotides6,MelissaT6],[MelissaNegativeHStickyNucleotides4,MelissaT4],[TakinoueNegativeHStickyNucleotides8,TakinoueT8],[TakinoueNegativeHStickyNucleotides,TakinoueT])


    """----------------------------ΔS----------------------------"""
    MelissaSequence6S = -118.0
    MelissaSequence4S = -76.0
    TakinoueSequence8S = -162.4
    TakinoueSequencesS = [-76.0,-113.3,-162.4,-201.9,-253.8]

    MelissaSStickyNucleotides6 = [-472.0,-708.0,-944.0,-1416.0,-2832.0]
    MelissaSStickyNucleotides4 = -304.0
    TakinoueNegativeSStickyNucleotides8 = [-487.2,-649.6,-974.4]
    TakinoueNegativeSStickyNucleotides = [-228.0,-339.9,-487.2,-605.7,-761.4]

    """Estimated transition temperatures for corresponding ΔS multiplied by # of arms"""
    #SStickyNucleotidesvsT([MelissaSStickyNucleotides6,MelissaT6],[MelissaSStickyNucleotides4,MelissaT4],[TakinoueNegativeSStickyNucleotides8,TakinoueT8],[TakinoueNegativeSStickyNucleotides,TakinoueT])


    """----------------------------ΔG----------------------------"""
    MelissaSequence6G = [-6.824045000000002,-5.092245000000004,-4.7211450000000035,-3.978945,-3.6078450000000046]
    MelissaSequence4G = -5.559445000000003
    TakinoueSequence8G = [-5.156985000000004,-4.4845850000000045,-3.475985000000005]
    TakinoueSequencesG = [-3.9009350000000045,-4.347450000000002,-5.375515000000005,-5.868180000000001,-6.874725000000007]


    MelissaGStickyNucleotides6 = [-27.296180000000007, -30.553470000000022, -37.76916000000003, -47.74734, -86.58828000000011]
    MelissaGStickyNucleotides4 = -22.23778000000001
    TakinoueGStickyNucleotides8 = [-15.470955000000012, -17.938340000000018, -20.85591000000003]
    TakinoueGStickyNucleotides = [-11.702805000000014, -13.042350000000006, -16.126545000000014, -17.60454, -20.624175000000022]

    GStickyNucleotidesvsT([MelissaGStickyNucleotides6,MelissaT6],[MelissaGStickyNucleotides4,MelissaT4],[TakinoueGStickyNucleotides8,TakinoueT8],[TakinoueGStickyNucleotides,TakinoueT])


def StickyNucleotidesvsT(data1,data2,data3,data4):
	plt.plot(data1[0],data1[1],label="Melissa's data (4,6,8,12,24 arms: 6 nt)",color='springgreen',marker='o',markersize=5,linestyle="None")
	plt.plot(data2[0],data2[1],label="Melissa's data (4 arms: 4 nt)",color="mediumturquoise",marker='o',markersize=5,linestyle="None")
	plt.plot(data3[0],data3[1],label="Takinoue's data (3,4,6 arms: 8 nt)",color="gold",marker='o',markersize=5,linestyle="None")
	plt.plot(data4[0],data4[1],label="Takinoue's data (3 arms: 4,6,8,10,12 nt)",color="orange",marker='o',markersize=5,linestyle="None")

	plt.xlabel("# of 'sticky nucleotides'")
	plt.ylabel("Temperature (°C)")
	plt.title("Estimated transition temperatures for corresponding # of 'sticky nucleotides'")
	plt.legend()

	plt.show()

def HStickyNucleotidesvsT(data1,data2,data3,data4):
	plt.plot(data1[0],data1[1],label="Melissa's data (4,6,8,12,24 arms: 6 nt)",color='springgreen',marker='o',markersize=5,linestyle="None")
	plt.plot(data2[0],data2[1],label="Melissa's data (4 arms: 4 nt)",color="mediumturquoise",marker='o',markersize=5,linestyle="None")
	plt.plot(data3[0],data3[1],label="Takinoue's data (3,4,6 arms: 8 nt)",color="gold",marker='o',markersize=5,linestyle="None")
	plt.plot(data4[0],data4[1],label="Takinoue's data (3 arms: 4,6,8,10,12 nt)",color="orange",marker='o',markersize=5,linestyle="None")

	plt.xlabel("ΔH (kcal/mol)")
	plt.ylabel("Temperature (°C)")
	plt.title("Estimated transition temperatures for corresponding ΔH multiplied by # of arms")
	plt.legend()

	plt.show()

def NegativeHStickyNucleotidesvsT(data1,data2,data3,data4):
    plt.plot(data1[0],data1[1],label="Melissa's data (4,6,8,12,24 arms: 6 nt)",color='springgreen',marker='o',markersize=5,linestyle="None")
    plt.plot(data2[0],data2[1],label="Melissa's data (4 arms: 4 nt)",color="mediumturquoise",marker='o',markersize=5,linestyle="None")
    plt.plot(data3[0],data3[1],label="Takinoue's data (3,4,6 arms: 8 nt)",color="gold",marker='o',markersize=5,linestyle="None")
    plt.plot(data4[0],data4[1],label="Takinoue's data (3 arms: 4,6,8,10,12 nt)",color="orange",marker='o',markersize=5,linestyle="None")

    plt.xlabel("-ΔH (kcal/mol)")
    plt.ylabel("Temperature (°C)")
    plt.title("Estimated transition temperatures for corresponding ΔH multiplied by # of arms")
    plt.legend()

    plt.show()

def SStickyNucleotidesvsT(data1,data2,data3,data4):
    plt.plot(data1[0],data1[1],label="Melissa's data (4,6,8,12,24 arms: 6 nt)",color='springgreen',marker='o',markersize=5,linestyle="None")
    plt.plot(data2[0],data2[1],label="Melissa's data (4 arms: 4 nt)",color="mediumturquoise",marker='o',markersize=5,linestyle="None")
    plt.plot(data3[0],data3[1],label="Takinoue's data (3,4,6 arms: 8 nt)",color="gold",marker='o',markersize=5,linestyle="None")
    plt.plot(data4[0],data4[1],label="Takinoue's data (3 arms: 4,6,8,10,12 nt)",color="orange",marker='o',markersize=5,linestyle="None")

    plt.xlabel("ΔS (cal mol-1 K-1)")
    plt.ylabel("Temperature (°C)")
    plt.title("Estimated transition temperatures for corresponding ΔS multiplied by # of arms")
    plt.legend()

    plt.show()

def GStickyNucleotidesvsT(data1,data2,data3,data4):
    plt.plot(data1[0],data1[1],label="Melissa's data (4,6,8,12,24 arms: 6 nt)",color='springgreen',marker='o',markersize=5,linestyle="None")
    plt.plot(data2[0],data2[1],label="Melissa's data (4 arms: 4 nt)",color="mediumturquoise",marker='o',markersize=5,linestyle="None")
    plt.plot(data3[0],data3[1],label="Takinoue's data (3,4,6 arms: 8 nt)",color="gold",marker='o',markersize=5,linestyle="None")
    plt.plot(data4[0],data4[1],label="Takinoue's data (3 arms: 4,6,8,10,12 nt)",color="orange",marker='o',markersize=5,linestyle="None")

    plt.xlabel("ΔG (kcal/mol)")
    plt.ylabel("Temperature (°C)")
    plt.title("Estimated transition temperatures for corresponding ΔG multiplied by # of arms")
    plt.legend()

    plt.show()

main()
