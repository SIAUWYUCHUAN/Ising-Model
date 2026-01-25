import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from Lattice_Class import Lattice, Observables
import numpy as np
import argparse

# python Animate.py -N 50 -J 1 -T 2 -iter 1000 -algo glauber -therm 10 -sampling 10 -config hot

parser = argparse.ArgumentParser(description="2D Ising Model Simulation")
parser.add_argument("-N", type=int, default=50, help="Lattice size N")
parser.add_argument("-J", type=float, default=1.0, help="Coupling constant")
parser.add_argument("-T", type=float, default=0.1, help="Temperature in units of J/kB")
parser.add_argument("-iter", type=int, default=1000, help="Number of iterations")
parser.add_argument("-algo", type=str, default="glauber", help="Algorithm: 'kawasaki' or 'glauber'")
parser.add_argument("-therm", type=int, default=10, help="Thermalisation steps")
parser.add_argument("-sampling", type=int, default=10, help="Sampling interval")
parser.add_argument("-config", type=str, default="hot", help="'hot' or 'cold' start")
parser.add_argument("-error", type=str, default="jackknife", help="Error analysis method: 'jackknife' or 'bootstrap'")

args = parser.parse_args()
ising = Lattice(size=args.N, J=args.J, T=args.T, iterations=args.iter,
                algorithm=args.algo, thermalisation=args.therm,
                sampling=args.sampling, start_config=args.config)

ani = ising.animate()

compute = input("Compute Observables after Thermalisation? (yes/no): ")

if compute == 'yes':
    Measurements = Observables(ising)
    Measurements.observables_with_errors(args.error)
    print("Average Magnetisation: ", Measurements.avg_mag, "±", Measurements.mag_error)
    print("Average Energy: ", Measurements.avg_E, "±", Measurements.E_error)  
    print("Susceptibility: ", Measurements.susceptibility, "±", Measurements.susceptibility_error)
    print("Heat Capacity: ", Measurements.heat_capacity, "±", Measurements.heat_capacity_error)


