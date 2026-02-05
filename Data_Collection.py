import matplotlib.pyplot as plt
from Lattice_Class_New import Lattice, Observables
import numpy as np
import json
import os
import argparse

# Uncomment the last part to run on terminal
# "python Data_Collection.py -size 15 -iter 1000 -sampling 10 -therm 200 -algo glauber -error jackknife -T_max 3.0 -T_min 1.0 -T_step 0.1 -direction heating -J 1 -plotting 1"

def iterate_T(output_dir, size, iteration, sampling, thermalisation, algorithm, error_analysis, T, J, old_grid):
    """
    Iterate T runs the ising model for a single temperature by first creating a ising object using the Lattice Class
    !!! The starting config here is a pre specified array from a previous simulation !!!

    """

    ising = Lattice(size=size, J=J, T=T, iterations=iteration,
                    algorithm=algorithm, thermalisation=thermalisation,
                    sampling=sampling, start_config=old_grid)
    
    list(ising.sim())  # Run the simulation
    new_grid = ising.grid

    # Compute observables
    Measurements = Observables(ising)
    Measurements.observables_with_errors(function=error_analysis)

    data = {
        "T": float(T),
        "magnetisation": [float(x) for x in ising.magnetisation],
        "total_energy": [float(x) for x in ising.totenergy]
    }
    filename = os.path.join(output_dir, f"T_{T:.3f}.json")

    with open(filename, "w") as f:
        json.dump(data, f, indent=4)

    return new_grid, Measurements


def plot(T_vals, mag_list, mag_list_err, E_list, E_list_err, susceptibility_list, susceptibility_list_err, heat_capacity_list, heat_capacity_list_err, size, algorithm):

    plt.figure(figsize=(11, 9))
    plt.suptitle(f"2D Ising Model Simulations for N = {size} under {algorithm} dynamics", fontsize=15)

    # Average Magnetization per spin
    plt.subplot(2, 2, 1)
    plt.errorbar(T_vals, mag_list, yerr=mag_list_err, fmt='o', markersize=2, elinewidth=2, capsize=3)
    plt.title(r"Average Magnetization $\langle |M| \rangle$ against Temperature [$k_bT$]", fontsize=10)
    plt.xlabel(r"Temperature [$k_bT$]", fontsize=9)
    plt.ylabel(r"Average Magnetization $\langle |M| \rangle$", fontsize=9)

    # Average Energy per spin
    plt.subplot(2, 2, 2)
    plt.errorbar(T_vals, E_list, yerr=E_list_err, fmt='s', markersize=2,
                elinewidth=2, capsize=3)
    plt.title(r"Average Energy  $\langle E \rangle$ [J] against Temperature [$k_bT$]", fontsize=10)
    plt.xlabel(r"Temperature [$k_bT$]", fontsize=9)
    plt.ylabel(r"Average Energy  $\langle E \rangle$ [J] ", fontsize=9)

    # Susceptibility
    plt.subplot(2, 2, 3)
    plt.errorbar(T_vals, susceptibility_list, yerr=susceptibility_list_err, fmt='^', markersize=2,
                elinewidth=2, capsize=3)
    plt.title(r"Susceptibility $\chi$ against Temperature [$k_bT$]", fontsize=10)
    plt.xlabel(r"Temperature [$k_bT$]", fontsize=9)
    plt.ylabel(r"Susceptibility $\chi$", fontsize=9)

    # Heat Capacity
    plt.subplot(2, 2, 4)
    plt.errorbar(T_vals, heat_capacity_list, yerr=heat_capacity_list_err, fmt='d', markersize=2,
                elinewidth=2, capsize=3)
    plt.title(r"Heat Capacity $C_v$ per spin [$k_b$] against Temperature [$k_bT$]", fontsize=10)
    plt.xlabel(r"Temperature [$k_bT$]", fontsize=9)
    plt.ylabel(r"Heat Capacity $C_v$ per spin [$k_b$]", fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()


def run_all(output_dir, size, iteration, sampling, thermalisation, algorithm, error_analysis, T_vals, J, grid, plotting):
    """
    Run_all computes the observables for Average magnetisation, Average energy, Heat capacity and susceptibility for a set of temperature values T_vals
    """
    mag_list, mag_list_err = [], []
    E_list, E_list_err = [], []
    susceptibility_list, susceptibility_list_err = [], []
    heat_capacity_list, heat_capacity_list_err = [], []

    for i, T in enumerate(T_vals):
        if i == 0:
            grid, Measurements = iterate_T(output_dir, size, iteration, sampling, thermalisation * 5, algorithm, error_analysis, T, J, grid)
        
        else:
            grid, Measurements = iterate_T(output_dir, size, iteration, sampling, thermalisation, algorithm, error_analysis, T, J, grid)

        mag_list.append(Measurements.avg_mag)
        mag_list_err.append(Measurements.mag_error)
        E_list.append(Measurements.avg_E)
        E_list_err.append(Measurements.E_error)
        susceptibility_list.append(Measurements.susceptibility)
        susceptibility_list_err.append(Measurements.susceptibility_error)
        heat_capacity_list.append(Measurements.heat_capacity)
        heat_capacity_list_err.append(Measurements.heat_capacity_error)

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
    
    if plotting == 1:
        plot(T_vals, mag_list, mag_list_err, E_list, E_list_err, susceptibility_list, susceptibility_list_err, heat_capacity_list, heat_capacity_list_err, size, algorithm)

    filename = os.path.join(output_dir, "Observables.json")
    with open(filename, "w") as f:
        json.dump(results, f, indent=4)

    print("All_Results are saved into", filename)


def collect_data(size, iteration, sampling, thermalisation, algorithm, error_analysis, direction, T_max, T_min, T_step, J, plotting):
    """
    Final funtion to initiate all data collection
    users can specify the direct of heat transfer (Cooling or heating) - which will results in different starting conditions
    """

    output_dir = f"{algorithm}_{size}_{direction}" # Output directory
    os.makedirs(output_dir, exist_ok=True)

    if direction == 'cooling':
        T_vals = np.arange(T_max, T_min - T_step, -T_step)

    elif direction == 'heating':
        T_vals = np.arange(T_min, T_max + T_step, T_step)

    if algorithm == 'kawasaki':
        grid = np.random.choice([-1, 1], size=(size, size)) # For kawasaki, always start with a random configuration

    elif algorithm == 'glauber':
        if direction == 'cooling':
            grid = np.random.choice([-1, 1], size=(size, size)) # random initial configuration for cooling
        elif direction == 'heating':
            grid = np.ones((size, size), dtype=int) # ordered spin configuration up for heating

    run_all(output_dir, size, iteration, sampling, thermalisation, algorithm, error_analysis, T_vals, J, grid, plotting)


#______________________________________________________________________________________________________________________________

## Comment this out to run it on terminal (Running on Jypter notebook is currently faster)
# -----------------------------
# Argument parser
# -----------------------------
# parser = argparse.ArgumentParser(description="2D Ising Model Simulation over Temperature Range")

# parser.add_argument("-size", type=int, default=50, help="Lattice size N")
# parser.add_argument("-iter", type=int, default=10000, help="Number of iterations per temperature")
# parser.add_argument("-sampling", type=int, default=10, help="Sampling interval")
# parser.add_argument("-therm", type=int, default=1000, help="Thermalisation steps")
# parser.add_argument("-algo", type=str, default='glauber', choices=['glauber','kawasaki'], help="Update algorithm")
# parser.add_argument("-error", type=str, default='jackknife', choices=['jackknife','bootstrap'], help="Error analysis method")
# parser.add_argument("-T_max", type=float, default=3.0, help="max temperature")
# parser.add_argument("-T_min", type=float, default=1.0, help="min temperature")
# parser.add_argument("-T_step", type=float, default=0.1, help="Temperature step_size")
# parser.add_argument("-direction", type=str, default='cooling', choices=['cooling','heating'], help="Direction of temperature change")
# parser.add_argument("-J", type=float, default= 1, help="Set the J value")
# parser.add_argument("-plotting", type=int, default= 1, help="key 1 for plotting and 0 for no plots")

# args = parser.parse_args()

# # -----------------------------
# # Set parameters
# # -----------------------------
# size = args.size
# iteration = args.iter
# sampling = args.sampling
# thermalisation = args.therm
# algorithm = args.algo
# error_analysis = args.error
# direction = args.direction
# T_max = args.T_max
# T_min = args.T_min
# T_step = args.T_step
# J = args.J
# plotting = args.plotting

# collect_data(size, iteration, sampling, thermalisation, algorithm, error_analysis, direction, T_max, T_min, T_step, J, plotting)