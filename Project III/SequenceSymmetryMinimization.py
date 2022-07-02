import math

"""This Python program checks the viability of
a designed DNA strand and its base pairing interactions
using Sequence Symmetry Minimization,
which is extremely useful for labs that artificially
design double-stranded DNA using software such as NUPACK"""

def main():
    print("Input your sequences by strand with the last bases of the sequence as the sticky ends.")
    number = int(input("Number of strands: "))
    stickyendnumber = int(input("Number of sticky ends: "))
    sequences = []
    for i in range(number):
        sequence = str(input("Sequence "+str(i+1)+": ")).upper()
        sequences.append(sequence)

    viable = checkATGC(sequences)


    """if viable:
        viable = checkSequenceSymmetry(sequences, stickyendnumber)"""

    if viable:
        print("Sequences are viable")
    else:
        print("Sequences are not viable")


def checkATGC(sequences):
    for i in sequences:
        for j in i:
            if (j != "A" and j != "T" and j != "G" and j != "C"):
                return False
    return True

def checkSequenceSymmetry(sequences, stickyendnumber):
    sequenceswithoutstickyend = []
    for i in sequences:
        sequence = i[0:len(i)-stickyendnumber]
        sequenceswithoutstickyend.append(sequence)

    complementarysequences = []
    for i in sequences:
        sequence = ""
        for j in i:
            if (j == "A"):
                sequence += "T"
            elif (j == "T"):
                sequence += "A"
            elif (j == "G"):
                sequence += "C"
            elif (j == "C"):
                sequence += "G"
        complementarysequences.append(sequence)




main()
