This folder consists of three python script to simulate the 2D ising model

Lattice_Class.py consists of two classes
    1. Lattice Class initialises a lattice with the following parameters:
        size (Int): Size of the lattice (size x size)
        J (float): Interaction strength
        T (float): Temperature
        iterations (Int): Number of iterations
        algorithm (str): 'glauber' or 'kawasaki'
        thermalisation (Int): Number of thermalisation sweeps
        sampling (Int): Sampling frequency
        start_config: 'hot', 'cold', or a numpy array

        During the simulation the following Observables are calculated:
            Total Magnetisation of the lattice
            Total Energy of the Lattice

    2. Observable class hold attributes and methods that compute the errors and oberservables of interest
        It inherits the Lattice class with the arrays of Magnetisation and Eneergy which thenm computes the susceptibioity and heat capacity
        Users can choose between the Jackknife or the bootstrap error


Animate.py produces animations and tracks the Absolsute magnetisation and total energy of the lattice dynamically. 
    Simulation grids are updated at every sweep.
    
    To run the simulation use the following prompt into the command line
    "python Animate.py -N 50 -J 1 -T 2 -iter 1000 -algo glauber -therm 10 -sampling 10 -config hot"

    usage: Animate.py [-h] [-N N] [-J J] [-T T] [-iter ITER] [-algo ALGO] [-therm THERM] [-sampling SAMPLING] [-config CONFIG]
                  [-error ERROR]

    2D Ising Model Simulation

    options:
    -h, --help          show this help message and exit
    -N N                Lattice size N
    -J J                Coupling constant
    -T T                Temperature in units of J/kB
    -iter ITER          Number of iterations
    -algo ALGO          Algorithm: 'kawasaki' or 'glauber'
    -therm THERM        Thermalisation steps
    -sampling SAMPLING  Sampling interval
    -config CONFIG      'hot' or 'cold' start
    -error ERROR        Error analysis method: 'jackknife' or 'bootstrap'


Data_Collection.py runs a set of temperatures and collects the data of all obersvables into a folder names {Algorithm}_{Lattice size}_{heating / cooling} with the following json files

    1. T_{Temperature}.json: it has three labels with 
        a. T - Temperature
        b. magnetisation - list of all the magnetisation values of all sampled grids
        c. total_energy - list of total energy values of all sampled grids
    2. Observables.json: It holds all the computed obersvables and their errors for each temperature value:
        a. "T_vals": T_vals.tolist(),
        b. "parameters": {
                "size": size,
                "J": 1,
                "iterations": iteration,
                "algorithm": algorithm,
                "thermalisation": thermalisation,
                "sampling": sampling,
                "error_analysis": error_analysis
            },
        c. "magnetization": {"values": mag_list, "errors": mag_list_err},
        d. "energy": {"values": E_list, "errors": E_list_err},
        e. "susceptibility": {"values": susceptibility_list, "errors": susceptibility_list_err},
        f. "heat_capacity": {"values": heat_capacity_list, "errors": heat_capacity_list_err}



    usage: Data_Collection.py [-h] [-size SIZE] [-iter ITER] [-sampling SAMPLING] [-therm THERM] [-algo {glauber,kawasaki}]
                          [-error {jackknife,bootstrap}] [-T_max T_MAX] [-T_min T_MIN] [-T_step T_STEP]
                          [-direction {cooling,heating}]

    To run the simulation use the following prompt into the command line
    "python Data_Collection.py -size 15 -iteration 1000 -sampling 10 -thermalisation 200 -algorithm kawasaki -error_analysis jackknife -T_max 3.0 -T_min 1.0 -T_step 0.1 -direction heating"

    2D Ising Model Simulation over Temperature Range

    options:
    -h, --help            show this help message and exit
    -size SIZE            Lattice size N
    -iter ITER            Number of iterations per temperature
    -sampling SAMPLING    Sampling interval
    -therm THERM          Thermalisation steps
    -algo {glauber,kawasaki}
                            Update algorithm
    -error {jackknife,bootstrap}
                            Error analysis method
    -T_max T_MAX          max temperature
    -T_min T_MIN          min temperature
    -T_step T_STEP        Temperature step_size
    -direction {cooling,heating}
                            Direction of temperature change


