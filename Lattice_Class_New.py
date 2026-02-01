import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Lattice:
    """ 
    Initialisation of the lattice requires the following parameters:
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
    """

    def __init__(self, size, J, T, iterations, algorithm='glauber', thermalisation = 100, sampling = 10, start_config=None):
        self.size = size

        if isinstance(start_config, np.ndarray):
            self.grid = start_config.copy()

        elif start_config == "cold":
            self.grid = -np.ones((size, size), dtype=int)

        elif start_config == "hot":
            self.grid = np.random.choice([-1, 1], size=(size, size))

        else:
            # default: random (hot start)
            self.grid = np.random.choice([-1, 1], size=(size, size))

        self.J = J
        self.T = T
        self.thermalisation = thermalisation
        self.sampling = sampling
        self.iteration = iterations
        self.algorithm = algorithm

        self.magnetisation = []
        self.totenergy = []

        # Compute the initial Energy and Magnetisation
        self.current_magnetisation = np.sum(self.grid) # to keep track of the current magnetisation

        total_energy = 0
        for i in range(self.size):
            for j in range(self.size):
                total_energy += self.local_energy(i, j)

        self.current_totenergy = total_energy / 2  # to keep track of the current total energy

        num_sites = (self.iteration + self.thermalisation) * self.size ** 2 # number of iterations

        # Generate all random numbers at one shot
        if self.algorithm == "glauber":
            self.i1_coords = np.random.randint(0, self.size, size=num_sites)
            self.j1_coords = np.random.randint(0, self.size, size=num_sites)
            self.k = np.random.rand(num_sites)

        elif self.algorithm == "kawasaki":
            self.i1_coords = np.random.randint(0, self.size, size=num_sites)
            self.j1_coords = np.random.randint(0, self.size, size=num_sites)
            self.i2_coords = np.random.randint(0, self.size, size=num_sites)
            self.j2_coords = np.random.randint(0, self.size, size=num_sites)
            self.k = np.random.rand(num_sites)


    def local_energy(self, i , j):
        """
        Calculate the local energy of the site at (i, j) using the 4 nearest neighbours with periodic boundary conditions
        E = - J * s_ij * (s_left + s_right + s_up + s_down)
        """
        left = self.grid[i, (j - 1) % self.size]
        right = self.grid[i, (j + 1) % self.size]
        up = self.grid[(i - 1) % self.size, j]
        down = self.grid[(i + 1) % self.size, j]
        return -self.J * self.grid[i, j] * (left + right + up + down)
    
    def glauber_update(self, i, j, k):
        """
        Perform a single Glauber update using the metropolis algorithm
            1. Calculate the change in energy ΔE if the spin at (i, j) is flipped using ΔE = -2 * local_energy
            2. If ΔE <= 0, accept the new configuration
            3. If ΔE > 0, accept the new configuration with probability exp(-ΔE / T), otherwise revert the spin
        """

        old_energy = self.local_energy(i, j)
        proposed_spin = -self.grid[i, j]
        self.grid[i, j] = proposed_spin
        delta_E = -2 * old_energy # Change in energy
        delta_M = 2 * proposed_spin  # Change in magnetisation

        if delta_E <= 0:
            # Accept the new configuration
            self.current_totenergy += delta_E
            self.current_magnetisation += delta_M
            return
        else:
            acceptance_prob = np.exp(-delta_E / self.T)
            if k < acceptance_prob:
                # Accept the new configuration
                self.current_totenergy += delta_E
                self.current_magnetisation += delta_M
                return
            else:
                # Reject the new configuration, revert the spin
                self.grid[i, j] = -proposed_spin

    def periodic_dist(self, a, b, L):
        return min(abs(a - b), L - abs(a - b))

    def kawasaki_update(self, i1, j1, i2, j2, k):
        """
        Perform a single Kawasaki update using the metropolis algorithm
            1. Calculate the change in energy ΔE if the spins at (i1, j1) and (i2, j2) are swapped
                a. CHECK 1: If they are the same spin, do nothing 
                b. CHECK 2: If they are adjacent, subtract the interaction energy between them 
                ΔE = -2 * (local_energy_1 + local_energy_2) + 4 * J
                c. If they are not adjacent: ΔE = -2 * (local_energy_1 + local_energy_2)
            2. If ΔE <= 0, accept the new configuration
            3. If ΔE > 0, accept the new configuration with probability exp(-ΔE / T), otherwise revert the spins
        """

        if self.grid[i1, j1] == self.grid[i2, j2]:
            return  # No point in swapping identical spins

        old_energy_1 = self.local_energy(i1, j1)
        old_energy_2 = self.local_energy(i2, j2)

        # Swap spins
        proposed_spin_1 = self.grid[i2, j2]
        proposed_spin_2 = self.grid[i1, j1]
        self.grid[i1, j1] = proposed_spin_1
        self.grid[i2, j2] = proposed_spin_2

        delta_E = -2 * (old_energy_1 + old_energy_2)

        dx = self.periodic_dist(i1, i2, self.size)
        dy = self.periodic_dist(j1, j2, self.size)
        is_neighbor = (dx + dy == 1)

        if is_neighbor:
             delta_E += 4 * self.J

        if delta_E <= 0:
            self.current_totenergy += delta_E
            return
        
        else:
            acceptance_prob = np.exp(-delta_E / self.T)
            if k < acceptance_prob:
                self.current_totenergy += delta_E
                return
            else:
                self.grid[i1, j1] = proposed_spin_2
                self.grid[i2, j2] = proposed_spin_1
        

    def sweep(self, iter):
        """
        Performs a single sweep of the Metropolis algorithm over the entire lattice (N^2 proposed updates)
        *** Random numbers are generated in a single array to reduce the computation time. (2 sets of random numbers )
        """
        updates = self.size ** 2

        if self.algorithm == 'glauber':
            i_coords = self.i1_coords[iter * updates : (iter + 1) * updates]
            j_coords = self.j1_coords[iter * updates : (iter + 1) * updates]
            k = self.k[iter * updates : (iter + 1) * updates]

            for new_i, new_j, k_val in zip(i_coords, j_coords, k):
                self.glauber_update(new_i, new_j, k_val)

        elif self.algorithm == 'kawasaki':
            i1_coords = self.i1_coords[iter * updates : (iter + 1) * updates]
            j1_coords = self.j1_coords[iter * updates : (iter + 1) * updates]
            i2_coords = self.i2_coords[iter * updates : (iter + 1) * updates]
            j2_coords = self.j2_coords[iter * updates : (iter + 1) * updates]
            k = self.k[iter * updates : (iter + 1) * updates]

            for i1, j1, i2, j2, k_val in zip(i1_coords, j1_coords, i2_coords, j2_coords, k):
                while i1 == i2 and j1 == j2: # Ensure distinct sites, if the sites are the same, resample
                    i2 = np.random.randint(0, self.size)
                    j2 = np.random.randint(0, self.size)
                self.kawasaki_update(i1, j1, i2, j2, k_val)

    def sim(self):
        """
        Runs a complete simulation of the Ising model
        Yields the current grid, magnetisation and total energy at each sampling step
        1. Perform thermalisation sweeps
        2. Perform iteration sweeps, sampling every 'sampling' sweeps
        3. Yield the sampled grid, magnetisation and total energy
        """
        for i in range(self.thermalisation):
            self.sweep(i)

        for i in range(self.iteration):
            self.sweep(i + self.thermalisation)
            if i % self.sampling == 0:
                self.magnetisation.append(self.current_magnetisation)
                self.totenergy.append(self.current_totenergy)
                yield self.grid.copy(), self.magnetisation.copy(), self.totenergy.copy()


    def animate(self, interval=100):
        """
        Animates the Ising model simulation and keeps track of the Magnetisation and Total Energy
        """
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))
        im = ax1.imshow(self.grid, cmap="coolwarm", vmin=-1, vmax=1, animated=True)
        ax1.set_title("Lattice")
        cbar = plt.colorbar(im, ax=ax1, ticks=[-1, 1])
        cbar.ax.set_yticklabels(['Spin Down', 'Spin Up'])
        
        line, = ax2.plot([], [], 'o-', markersize=2)
        ax2.set_xlim(0, self.iteration)
        ax2.set_ylim(-5000, 5000)
        ax2.set_title("Magnetization")
        ax2.set_xlabel("Sweep")
        ax2.set_ylabel("M")

        line3, = ax3.plot([], [], 'o-', markersize=2)
        ax3.set_xlim(0, self.iteration)
        ax3.set_ylim(-5000, 0)
        ax3.set_title("Total Energy")
        ax3.set_xlabel("Sweep")
        ax3.set_ylabel("E")

        def update(frame):
            grid, mag, E = frame
            im.set_data(grid)
            x_vals = [k * self.sampling for k in range(len(mag))]
            line.set_data(x_vals, mag)
            line3.set_data(x_vals, E)
            return [im, line, line3]
        
        ani = FuncAnimation(fig, update, frames=self.sim(), interval=interval, blit=True)
        plt.tight_layout()
        plt.show()
        return ani
    
class Observables:
    def __init__(self, Lattice):
        """
        Initialisation of the Observables class requires a Lattice object
        It inherits the magnetisation and total energy data from the Lattice object
        """
        self.Lattice = Lattice
        self.magnetisation = Lattice.magnetisation
        self.totenergy = Lattice.totenergy

    def com_average_magnetisation(self, mag):
        return np.mean(np.abs(mag)) # it computes the absolute value of the magnetisation
    
    def com_average_energy(self, E):
        return np.mean(E)
    
    def com_susceptibility(self, mag):
        return (np.mean(np.array(mag) ** 2) - np.mean(mag) ** 2) / (self.Lattice.size**2 * self.Lattice.T)
    
    def com_heat_capacity(self, E):
        return (np.mean(np.array(E) ** 2) - np.mean(E) ** 2) / (self.Lattice.size**2 * self.Lattice.T ** 2)

    def observables(self):
        """
        Computes the average Magnetisation, average Energy, Susceptibility and Heat Capacity for an entire set of sampled states
        """
        self.susceptibility = self.com_susceptibility(self.magnetisation)
        self.heat_capacity = self.com_heat_capacity(self.totenergy)
        self.avg_mag = self.com_average_magnetisation(self.magnetisation)
        self.avg_E = self.com_average_energy(self.totenergy)

        print("Average Magnetization:", self.avg_mag, "Average Energy:", self.avg_E, "Susceptibility:", self.susceptibility, "Heat Capacity:", self.heat_capacity)
        return self.avg_mag, self.avg_E, self.susceptibility, self.heat_capacity
    
    def Jackknife_errors(self, arr, function):
        """
        Computes the Jackknife error for a set of data 'arr' using the statistic defined by 'function'
        1. computes the function on the entire array
        2. computes the function on n jackknife samples, each with one element removed
        3. computes the jackknife error using the formula
           jackknife = sqrt( sum (c_i - c)^2 )
        """
        n = len(arr)
        c = function(arr)
        jackknife_samples = np.array([function(np.delete(arr, i)) for i in range(n)])
        jackknife_error = np.sqrt(np.sum((ci - c)**2 for ci in jackknife_samples))
        return jackknife_error
    
    def Bootstrap_errors(self, arr, function, num_samples=1000):
        """
        Computes the Bootstrap error for a set of data 'arr' using the statistic defined by 'function'
        1. generates 'num_samples' bootstrap samples by resampling with replacement
        2. computes the function on each bootstrap sample
        3. computes the bootstrap error as the standard deviation of the bootstrap statistics
        """
        n = len(arr)
        bootstrap_samples = np.random.choice(arr, size=(num_samples, n), replace=True)
        bootstrap_statistics = np.array([function(sample) for sample in bootstrap_samples])
        bootstrap_error = np.std(bootstrap_statistics)
        return bootstrap_error
    
    def observables_with_errors(self, function = 'jackknife'):
        self.observables()

        if function == 'jackknife':
            err_function = self.Jackknife_errors
        elif function == 'bootstrap':
            err_function = self.Bootstrap_errors
        
        self.mag_error = err_function(self.magnetisation, self.com_average_magnetisation)
        self.E_error = err_function(self.totenergy, self.com_average_energy)
        self.susceptibility_error = err_function(self.magnetisation, self.com_susceptibility)
        self.heat_capacity_error = err_function(self.totenergy, self.com_heat_capacity)

        return self.mag_error, self.E_error, self.susceptibility_error, self.heat_capacity_error