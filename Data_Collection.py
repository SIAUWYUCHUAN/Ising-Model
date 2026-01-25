import matplotlib.pyplot as plt
from Lattice_Class import Lattice, Observables
import numpy as np
import json
import os
import argparse

# -----------------------------
# Argument parser
# -----------------------------
parser = argparse.ArgumentParser(description="2D Ising Model Simulation over Temperature Range")

parser.add_argument("-size", type=int, default=50, help="Lattice size N")
parser.add_argument("-iter", type=int, default=10000, help="Number of iterations per temperature")
parser.add_argument("-sampling", type=int, default=10, help="Sampling interval")
parser.add_argument("-therm", type=int, default=1000, help="Thermalisation steps")
parser.add_argument("-algo", type=str, default='glauber', choices=['glauber','kawasaki'], help="Update algorithm")
parser.add_argument("-error", type=str, default='jackknife', choices=['jackknife','bootstrap'], help="Error analysis method")
parser.add_argument("-T_max", type=float, default=3.0, help="max temperature")
parser.add_argument("-T_min", type=float, default=1.0, help="min temperature")
parser.add_argument("-T_step", type=float, default=0.1, help="Temperature step_size")
parser.add_argument("-direction", type=str, default='cooling', choices=['cooling','heating'], help="Direction of temperature change")

args = parser.parse_args()

# -----------------------------
# Set parameters
# -----------------------------
size = args.size
iteration = args.iter
sampling = args.sampling
thermalisation = args.therm
algorithm = args.algo
error_analysis = args.error
direction = args.direction
T_max = args.T_max
T_min = args.T_min
T_step = args.T_step

# Output directory
output_dir = f"{algorithm}_{size}_{direction}"
os.makedirs(output_dir, exist_ok=True)

# Temperature values
if direction == 'cooling':
    T_vals = np.arange(T_max, T_min - T_step, -T_step)

elif direction == 'heating':
    T_vals = np.arange(T_min, T_max + T_step, T_step)

# -----------------------------
# Lists to store results
# -----------------------------
mag_list, mag_list_err = [], []
E_list, E_list_err = [], []
susceptibility_list, susceptibility_list_err = [], []
heat_capacity_list, heat_capacity_list_err = [], []

# -----------------------------
# Function to update lattice at each temperature
# -----------------------------
def update(T, old_grid, thermalisation):
    ising = Lattice(size=size, J=1, T=T, iterations=iteration,
                    algorithm=algorithm, thermalisation=thermalisation,
                    sampling=sampling, start_config=old_grid)
    
    list(ising.sim())  # Run the simulation
    new_grid = ising.grid

    # Compute observables
    Measurements = Observables(ising)
    Measurements.observables_with_errors(function=error_analysis)

    # Store results
    mag_list.append(Measurements.avg_mag)
    mag_list_err.append(Measurements.mag_error)
    E_list.append(Measurements.avg_E)
    E_list_err.append(Measurements.E_error)
    susceptibility_list.append(Measurements.susceptibility)
    susceptibility_list_err.append(Measurements.susceptibility_error)
    heat_capacity_list.append(Measurements.heat_capacity)
    heat_capacity_list_err.append(Measurements.heat_capacity_error)

    # Save individual temperature JSON
    data = {
        "T": float(T),
        "magnetisation": [float(x) for x in ising.magnetisation],
        "total_energy": [float(x) for x in ising.totenergy]
    }
    filename = os.path.join(output_dir, f"T_{T:.3f}.json")
    with open(filename, "w") as f:
        json.dump(data, f, indent=4)

    return new_grid

# -----------------------------
# Run simulation over temperatures and initialize with grid type
# -----------------------------
if algorithm == 'kawasaki':
    old_grid = np.random.choice([-1, 1], size=(size, size)) # For kawasaki, always start with a random configuration

elif algorithm == 'glauber':
    if direction == 'cooling':
        old_grid = np.random.choice([-1, 1], size=(size, size)) # random initial configuration for cooling
    elif direction == 'heating':
        old_grid = np.ones((size, size), dtype=int) # ordered spin configuration up for heating


old_grid = update(T_vals[0], old_grid, thermalisation * 5) # Extra thermalisation for first T
for T in T_vals[1:]:
    old_grid = update(T, old_grid, thermalisation)


# -----------------------------
# Save all observables in one JSON
# -----------------------------
results = {
    "T_vals": T_vals.tolist(),
    "parameters": {
        "size": size,
        "J": 1,
        "iterations": iteration,
        "algorithm": algorithm,
        "thermalisation": thermalisation,
        "sampling": sampling,
        "error_analysis": error_analysis
    },
    "magnetization": {"values": mag_list, "errors": mag_list_err},
    "energy": {"values": E_list, "errors": E_list_err},
    "susceptibility": {"values": susceptibility_list, "errors": susceptibility_list_err},
    "heat_capacity": {"values": heat_capacity_list, "errors": heat_capacity_list_err}
}

filename = os.path.join(output_dir, "Observables.json")
with open(filename, "w") as f:
    json.dump(results, f, indent=4)

print(f"Results saved to {filename}")
