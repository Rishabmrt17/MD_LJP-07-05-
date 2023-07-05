#!/usr/bin/env python
# coding: utf-8

# In[5517]:


#import modules
import numpy as np
import matplotlib.pyplot as plt
import math


# In[5518]:


def run_md(N_atoms, r, T, m, size, step, delta_step, dt):
    """
    MD program using velocity verlet algorithm
    N_atoms = number_of_atoms
    r = distance_between_particles
    T = Boltzmann_temperature
    m = mass_of_particle
    size = size_of_the_box
    step = number_of_steps
    delta_step = frequency_in_steps
    dt = time_step
    """
    # Main MD loop
    # open trajectory file
    file = open("traj_MD.xyz", 'w')
        # initialize positions
    positions = initialize_positions(N_atoms, size)
    
        # initialize velocities
    velocities = initialize_velocities(m, kB, T)
    
    for steps in range(step):
    
        # Propagate Positions
        update_positions(positions, velocities, f, dt, size)
        
        # Propagate Velocities
        update_velocities(velocities, f, dt)
        
        if (step%delta_step==0):
            write_trajectory_frame(positions, file, step)


# In[5519]:


#################Sub-Routines#################

# initialize positions
def initialize_positions(N_atoms, size):
    """Initialize positions"""
    positions = np.random.rand(N_atoms, 3) * size
    
    return positions


# In[5520]:


N_atoms = 2
size = 10
positions = np.random.rand(N_atoms, 3) * size
print(positions)


# In[5521]:


def lj_potential(r, epsilon, sigma):
    """Lennard-Jones potential function"""
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)


# In[5522]:


def compute_distance(positions):
    """Compute distance between two particles"""
    return np.linalg.norm(positions[1] - positions[0])


# In[5523]:


# Compute forces
forces = np.zeros_like(positions)
r = compute_distance(positions)
f = lj_potential(r, epsilon, sigma) / r
forces[0] += f
forces[1] -= f


# In[5524]:


print(f)


# In[ ]:





# In[5525]:


# initialize velocities
def initialize_velocities(m, kB, T):
    """
    kB = boltzmann_temp
    """
    velocities = np.random.normal(loc=0, scale=1, size=(N_atoms, 3)) * np.sqrt(kB*T/m)
    
    return velocities


# In[5526]:


N_atoms = 2
kB = 0.08314
T = 300
m = 12
velocities = np.random.normal(loc=0, scale=np.sqrt(T), size=(N_atoms,3))
print(velocities)


# In[5527]:


# Propagate Positions
def update_positions(positions, velocities, dt, f, size):
    """
    positions = particle_positions
    velocities = particle_velocities
    dt = time_step
    f = force_on_particles
    size = size_of_box
    """
    
    positions += velocities * dt + 0.5 * f * dt**2
    
     # wrap into central box (box is from 0 to size in each dimension)
    for i in range(N):
        for j in range(N):
            if positions[i,j] < 0:
                positions[i,j] += size
            elif positions[i,j] > size:
                positions[i,j] -= size
    return positions


# In[5528]:


dt = 0.01
size = 10
positions += velocities * dt + 0.5 * f * dt**2


# In[5529]:


print(positions)


# In[5530]:


# Propagate Velocities
def update_velocities(velocities, f, dt):
    """
    velocities = particle_velocities
    f = forces
    dt = time_step
    
    """
    
    velocities += f * dt
    return velocities


# In[5531]:


dt = 0.01
f = 24 * epsilon / sigma * (2 * (sigma / r) ** 13 - (sigma / r) ** 7)
velocities += f * dt
print(velocities)


# In[5532]:


# Compute Energy
def kinetic_energy(N_atoms, kB, T):
    """kB - Boltzmann constant"""
    
    Energy = 1.5*N_atoms*kB*T
    return Energy


# In[5533]:


N_atoms = 2
kB = 0.08314
T = 300
Energy = 1.5*N_atoms*kB*T
print(Energy)


# In[5544]:


# Trajectory frame        
def write_trajectory_frame(positions, file_pointer, step):
    """
    positions = particle_positions
    file_pointer = trajectory_file_pointer
    step = step_number
    """
    for i in range(N_atoms):
        file_pointer.write("C %10.5f %10.5f %10.5f\n" % ( positions[i,0],  positions[i,1], positions[i,2]))
    
    # close trajectory file
    file.close()


# In[5545]:


#test sim
N_atoms = 2
r = sigma
T = 300
m = 12
size = 10
step = 100
delta_step = 10
dt = 0.01
sim = run_md(N_atoms, r, T, m, size, step, delta_step, dt)


# In[5546]:


print(sim)


# In[5547]:


print(file)


# In[5548]:


plt.plot(velocities)


# In[5549]:


plt.plot(positions)


# In[5561]:


def compute_gr(positions,size, N_atoms, dr, density):
    """
    Compute the radial distribution function, g(r), for a given set of particle positions.

    
    positions : Array of particle positions of shape (N_atoms, 3).
    size : Length of the simulation box.
    N_atoms : Number of particles in the system.
    dr : Width of the distance bins for g(r) calculation.
    density : Density of the system.
    """
    rdf = np.zeros(int(size / 2 / dr))
    volume = size ** 3
    N_atoms = len(positions)

    for i in range(N_atoms):
        for j in range(i + 1, N_atoms):
            rij = positions[i] - positions[j]
            rij -= size * np.round(rij / size)  # Apply periodic boundary conditions
            r = np.linalg.norm(rij)

            if r < size / 2:
                rdf_idx = int(r / dr)
                rdf[rdf_idx] += 1

    # Normalize g(r)
    num_pairs = N_atoms * (N_atoms - 1) / 2
    normalization_factor = 4 * np.pi * density * dr * num_pairs * N_atoms / volume
    rdf /= normalization_factor

    r = np.arange(0, size / 2, dr)
    return r, rdf


# In[5573]:


# Simulation parameters
size = 10.0
N_atoms = 2
dr = 0.1

# Varying densities
densities = [0.002, 0.5, 1.0]  # Example densities

# Perform simulations at different densities
for density in densities:
    # Compute number of particles based on density
    N_atoms = int(density * size ** 3)

    # Compute g(r) for the current density
    r, gr = compute_gr(positions, size, N_atoms, dr, density)


# In[5574]:


print(gr)


# In[5575]:


print(r)


# In[5576]:


plt.plot(r, gr)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[1090]:





# In[1091]:





# In[ ]:




