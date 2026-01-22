import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from Lattice_Class import Lattice
import numpy as np

size = 50
ising = Lattice(size=size, J=1, T=1, iterations=1000, algorithm='glauber', thermalisation = 10, sampling = 10)
ising.animate()
ising.observables()

ising = Lattice(size=size, J=1, T=1.5, iterations=1000, algorithm='kawasaki', thermalisation = 10, sampling = 10)
ising.animate()
ising.observables()

# ising = Lattice(size=size, J=1, T=2, iterations=1000, algorithm='glauber', thermalisation = 100, sampling = 10)
# ising.animate()
# ising.observables()

# ising = Lattice(size=size, J=1, T=2.5, iterations=1000, algorithm='glauber', thermalisation = 100, sampling = 10)
# ising.animate()
# ising.observables()

# ising = Lattice(size=size, J=1, T=3, iterations=1000, algorithm='glauber', thermalisation = 100, sampling = 10)
# ising.animate()
# ising.observables()