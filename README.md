# SummerResearch

During the summer of 2020, I worked as an Undergraduate Research Fellow at the Brandeis Materials Research Science and Engineering Center. In the Rogers Lab under PI Dr. Ben Rogers and postdoc Dr. Melissa Rinaldin, I conducted research on nanostructures called DNA nanostars with funding I acquired from the Summer Materials Undergraduate Research Fellowship provided by the National Science Foundation.

The following are projects I coded myself that aided the research:


## Project I
Project I explores the thermodynamic properties of DNA nanostars based on the nearest-neighbor (NN) model for nucleic acids.
- <b>StickyEndThermodynamics.py</b> takes in the sequence and transition temperature to output ∆G.
- <b>StickyEndThermodynamicsExtension.py</b> is a raw version of a large collection of graphs displaying the thermodynamic properties of these nanostars using data from Takinoue (2020) and my mentor Melissa's data. The graphs are currently being commented out, so need to remove '#' prior running the code. Also need to install numpy and matplotlib prior running the code.
- <b>NNModel.py</b> is the more user-friendly version of the first part of StickyEndThermodynamicsExtension.py, which involves the fundamental thermodynamic properties. The graphs are currently being commented out, so need to remove '#' prior running the code. Also need to install numpy and matplotlib prior running the code.
- <b>NNModel2.py</b> is the more user-friendly version of the second part of StickyEndThermodynamicsExtension.py, which involves accounting for 'sticky nucleotides' or the product of # of arms and # of sticky ends. The graphs are currently being commented out, so need to remove '#' prior running the code. Also need to install numpy and matplotlib prior running the code.

## Project II
Project II complements Project I by specifically investigating bond probability and particle valence using Wertheim’s theory. <b>WertheimTheory.py</b> takes in a list of bond probabilities and particle valences and outputs free energy estimations.

## Project III
Project III checks the viability of a designed DNA strand and its base pairing interactions using Sequence Symmetry Minimization, which is extremely useful for labs that artificially design double-stranded DNA using software such as NUPACK. <b>SequenceSymmetryMinimization.py</b> takes in the sequence and outputs its viability.
