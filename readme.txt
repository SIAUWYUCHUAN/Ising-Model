This folder consists of three python script to simulate the 2D ising model.

(A) Lattice_Class_new.py is a class file consisting of classes:
    1. Lattice Class initialises a lattice with the following initialisation attributes:
        size (Int): Size of the lattice (size x size)
        J (float): Interaction strength
        T (float): Temperature in kbT
        iterations (Int): Total number of iterations to during after equilibration
        algorithm (str): 'glauber' or 'kawasaki' 
        thermalisation (Int): Number of thermalisation sweeps for equilibration
        sampling (Int): Sampling frequency, number of iteration to be run before taking a energy or magnetisation measurement
        start_config: 'hot' - disordered array, 'cold' - regular array of all spin up, or a specified starting configuration in a numpy array

        General methods:
            It runs #(thermalisation + iteration) sweeps of a specified (Size, J, T) configuration for a specified algorithm (glauber or kawasaki)
            Each sweep attempt size**2 number of metropolis updates and samples the configuration for the magnetisation and energy at every interval of the sampling frequency

        During the simulation the following Observables are calculated and saved as attributes of the lattice class which is passed into the Observable class:
            Total Magnetisation of the lattice
            Total Energy of the Lattice


    2. Observable class hold attributes and methods that compute the errors and oberservables of interest
        It inherits the Lattice class with the arrays of Magnetisation and Energy which thenm computes the susceptibioity and heat capacity
        Users can choose between the Jackknife or the bootstrap error

_________________________________________________________________________________________________________________________

(B) Animate.py produces animations and tracks the Absolsute magnetisation and total energy of the lattice dynamically. 
    1. Simulation grids and data points for magnetisation and energy are updated at every sampled sweep.
    2. Simulation prompts user if they want to compute the full set of observables and their errors
            An example output is as follows:
            Compute Observables after Thermalisation? (yes/no): yes
            Average Magnetization: 2384.86 Average Energy: -4829.4 Susceptibility: 22.54761077333316 Heat Capacity: 4.524053333334128
            Average Magnetisation:  2384.86 ± 29.37180957368644
            Average Energy:  -4829.4 ± 16.113501611704297
            Susceptibility:  22.54761077333316 ± 10.117100578707145
            Heat Capacity:  4.524053333334128 ± 2.1110598254503703
            
    To run the simulation copy the following prompt into the command line
    -- "python Animate.py -N 50 -J 1 -T 2 -iter 1000 -algo glauber -therm 10 -sampling 10 -config hot" --

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

_________________________________________________________________________________________________________________________


(C) Data_Collection.py holds a set of function that are used to run the monte carlo simulation for a set of temperatures 

A single run of the collect_data() function outputs:
    - plots of the Average absolute magnetisation, average energy, heat capacity and the susceptibility against temperature [kbT]
    - folder names {Algorithm}_{Lattice size}_{heating / cooling} with the following json files
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


This file can be run in two ways:
    1. imported into a jypter notebook to be run (This is somehow faster)
            Open the jypter notebook Test_cases_jypter.ipynb and edit the parameters to the desired start parameters before running
            
            # TEST CASE
            size = 15
            iteration = 1000
            sampling = 10
            thermalisation = 100
            T_max = 3.0
            T_min = 1.0
            T_step = 0.1
            error_analysis = "bootstrap"  # "jackknife" or "bootstrap"
            J = 1
            plotting = 1

            algorithm = "glauber"  # "glauber" or "kawasaki"
            direction = "cooling"

            collect_data(size, iteration, sampling, thermalisation, algorithm, error_analysis, direction, T_max, T_min, T_step, J, plotting)

    2. To run on the terminal uncomment out the final section with the text (# Comment this out to run it on terminal)

        To run the simulation use the following prompt into the command line
        -- "python Data_Collection.py -size 15 -iter 1000 -sampling 10 -therm 200 -algo glauber -error jackknife -T_max 3.0 -T_min 1.0 -T_step 0.1 -direction heating -J 1 -plotting 1" --

        usage: Data_Collection.py [-h] [-size SIZE] [-iter ITER] [-sampling SAMPLING] [-therm THERM] [-algo {glauber,kawasaki}]
                            [-error {jackknife,bootstrap}] [-T_max T_MAX] [-T_min T_MIN] [-T_step T_STEP]
                            [-direction {cooling,heating}] [-J J] [-plotting PLOTTING]

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
        -J J                  Set the J value
        -plotting PLOTTING    key 1 for plotting and 0 for no plots

_________________________________________________________________________________________________________________________

(D) data_visualisation.ipynb plots the data from the kawasaki and glauber algorithm together to visually inspect the values
It also provides the values for all the raw data for each case

