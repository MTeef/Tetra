from tetralen import tetraLen
from numpy import array

# Tetrahedreon base lengths (Distances between Object-Points)
lab = 4
lbc = 3.2
lca = 3.1

# Head angles (Optical Center to Image Points angles)
ph1 = 0.5054
ph2 = 0.473
ph3 = 0.39
precision = .1

solutions, numberOfSolutions, elapsed = tetraLen(lab, lbc, lca, ph1, ph2, ph3, precision)
print (str(solutions) + "\n" + str(numberOfSolutions) + " Solutions found!")


