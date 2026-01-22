import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


class Lattice:
    def __init__(self, size, J, T, iterations, algorithm='glauber', thermalisation = 100, sampling = 10):
        self.size = size
        self.grid = np.random.choice([-1, 1], size=(size, size))
        self.J = J
        self.T = T
        self.thermalisation = thermalisation
        self.sampling = sampling
        self.iteration = iterations
        self.algorithm = algorithm
        self.magnetisation = []
        self.susceptibility = []
        self.totenergy = []

    def energy(self, i , j):
        # Calculate the energy of the site at (i, j)
        left = self.grid[i, (j - 1) % self.size]
        right = self.grid[i, (j + 1) % self.size]
        up = self.grid[(i - 1) % self.size, j]
        down = self.grid[(i + 1) % self.size, j]
        return -self.J * self.grid[i, j] * (left + right + up + down)
    
    def glauber_update(self, i, j):
        old_energy = self.energy(i, j)
        proposed_spin = -self.grid[i, j]
        self.grid[i, j] = proposed_spin
        new_energy = self.energy(i, j)
        delta_E = new_energy - old_energy

        if delta_E <= 0:
            # Accept the new configuration
            return
        else:
            acceptance_prob = np.exp(-delta_E / self.T)
            if rnd.random() < acceptance_prob:
                # Accept the new configuration
                return
            else:
                # Reject the new configuration, revert the spin
                self.grid[i, j] = -proposed_spin

    def kawaski_update(self, i1, j1, i2, j2):
        old_energy_1 = self.energy(i1, j1)
        old_energy_2 = self.energy(i2, j2)

        if self.grid[i1, j1] == self.grid[i2, j2]:
            return  # No point in swapping identical spins

        # Swap spins
        proposed_spin_1 = self.grid[i2, j2]
        proposed_spin_2 = self.grid[i1, j1]
        self.grid[i1, j1] = proposed_spin_1
        self.grid[i2, j2] = proposed_spin_2

        new_energy_1 = self.energy(i1, j1)
        new_energy_2 = self.energy(i2, j2)

        #check if sites are adjacent by checking the distance vector
        #why is it minus 4J, because we are removing the interaction between the two swapped spins
        if abs(i1 - i2) + abs(j1 - j2) == 1:
            delta_E = (new_energy_1 + new_energy_2) - (old_energy_1 + old_energy_2) - 4 * self.J * proposed_spin_1 * proposed_spin_2

        else:
            delta_E = (new_energy_1 + new_energy_2) - (old_energy_1 + old_energy_2)

        if delta_E <= 0:
            # Accept the new configuration
            return
        else:
            acceptance_prob = np.exp(-delta_E / self.T)
            if rnd.random() < acceptance_prob:
                # Accept the new configuration
                return
            else:
                # Reject the new configuration, revert the spins
                self.grid[i1, j1] = proposed_spin_2
                self.grid[i2, j2] = proposed_spin_1
        
    def sweep(self):
        num_sites = self.size ** 2

        if self.algorithm == 'glauber':
            # Generate all random i, j indices at once
            i_coords = np.random.randint(0, self.size, size=num_sites)
            j_coords = np.random.randint(0, self.size, size=num_sites)

            # Loop through the pre-generated coordinates
            for new_i, new_j in zip(i_coords, j_coords):
                self.glauber_update(new_i, new_j)

        elif self.algorithm == 'kawasaki':
            # Generate all random i1, j1, i2, j2 indices at once
            i1_coords = np.random.randint(0, self.size, size=num_sites)
            j1_coords = np.random.randint(0, self.size, size=num_sites)
            i2_coords = np.random.randint(0, self.size, size=num_sites)
            j2_coords = np.random.randint(0, self.size, size=num_sites)

            # Loop through the pre-generated coordinates
            for i1, j1, i2, j2 in zip(i1_coords, j1_coords, i2_coords, j2_coords):
                # Ensure distinct sites
                while i1 == i2 and j1 == j2:
                    i2 = np.random.randint(0, self.size)
                    j2 = np.random.randint(0, self.size)
                self.kawaski_update(i1, j1, i2, j2)


    def cal_magnetisation(self):
        self.magnetisation.append(np.sum(self.grid))

    def cal_total_energy(self):
        total_energy = 0
        for i in range(self.size):
            for j in range(self.size):
                total_energy += self.energy(i, j)
        self.totenergy.append(total_energy / 2)  # Each pair counted twice

    def observables(self):
        avg_mag = np.mean(np.abs(self.magnetisation))
        avg_E = np.mean(self.totenergy)

        susceptibility = (np.mean(np.array(self.magnetisation) ** 2) - np.mean(self.magnetisation) ** 2) / (self.size**2 * self.T)

        E_sq = np.mean(np.array(self.totenergy) ** 2)
        heat_capacity = (E_sq - avg_E ** 2) / (self.size**2 * self.T ** 2)
        print("Average Magnetization:", avg_mag, "Average Energy:", avg_E, "Susceptibility:", susceptibility, "Heat Capacity:", heat_capacity)

        return avg_mag, avg_E, susceptibility, heat_capacity

    def sim(self):
        self.magnetisation = []
        self.totenergy = []
        for i in range(self.thermalisation):
            self.sweep()
        
        for i in range(self.iteration):
            self.sweep()
            if i % self.sampling == 0:
                self.cal_magnetisation()
                self.cal_total_energy()
                yield self.grid.copy(), self.magnetisation.copy(), self.totenergy.copy()
        

    def animate(self, interval=100):
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))
        im = ax1.imshow(self.grid, cmap="coolwarm", vmin=-1, vmax=1, animated=True)
        ax1.set_title("Lattice")
        
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

