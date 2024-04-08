import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from .particles import Particles
from numba import jit, njit, prange, set_num_threads

"""
The N-Body Simulator class is responsible for simulating the motion of N bodies



"""

class NBodySimulator:

    def __init__(self, particles: Particles):
        
        # TODO
        self.particles = particles 
        self.setup() #default setup

        return

    def setup(self, G=1,
                    rsoft=0.01,
                    method="RK4",
                    io_freq=10,
                    io_header="nbody",
                    io_screen=True,
                    visualization=False):
        """
        Customize the simulation enviroments.

        :param G: the graivtational constant
        :param rsoft: float, a soften length
        :param meothd: string, the numerical scheme
                       support "Euler", "RK2", and "RK4"

        :param io_freq: int, the frequency to outupt data.
                        io_freq <=0 for no output. 
        :param io_header: the output header
        :param io_screen: print message on screen or not.
        :param visualization: on the fly visualization or not. 
        """
        
        # TODO
        self.G = G
        self.rsoft = rsoft
        self.io_freq = io_freq
        self.io_header = io_header
        self.io_screen = io_screen
        self.visualization = visualization
        if method.lower() == 'euler':
            _advance_particles = self._advance_particles_Euler()
        elif method.lower() == 'rk2':
            _advance_particles = self._advance_particles_RK2()
        elif method.lower() == 'rk4':
            _advance_particles = self._advance_particles_RK4()
        else:
            print('wrong method')
            quit() 



        return _advance_particles

    def evolve(self, dt:float, tmax:float):
        """
        Start to evolve the system

        :param dt: float, the time step
        :param tmax: float, the total time to evolve
        
        """

        # TODO
        time = self.particles.time
        dt = 0.01
        tmax = 50
        nstep = (tmax-time)/dt

        for n in range(int(nstep)):
            
            if (time + dt) > tmax:
                dt = tmax - time
            
            #update
            # check 10
            if (n % self.io_freq) == 0:

                fn = "data_"+self.io_header+"_"+str(n).zfill(6)+'.dat'
                fn = folder+"/"+fn







        print("Simulation is done!")
        return

    def _calculate_acceleration(self, nparticles, masses, positions):
        """
        Calculate the acceleration of the particles
        """
        accelerations = np.zeros_like(positions)
        nparticles = self.particles.nparticles

        # TODO
        for i in range(len(nparticles)):
            for j in range(nparticles):
                if j > i:
                    rij = positions[i]-positions[j]
                    r   = np.sqrt(np.sum(rij**2+self.rsoft**2))
                    force = -self.G*masses[i,0] * masses[j,0] /(r**3*rij) 




        return accelerations
        
    def _advance_particles_Euler(self, dt, particles):

        #TODO
        pos = particles.positions
        vel = particles.velocities
        acc = self._calculate_acceleration(self.particles.nparticles,self.particles.masses,pos)

        pos = pos+vel*dt
        vel = vel+acc*dt
        acc = self._calculate_acceleration(self.particles.nparticles,self.particles.masses,pos)

        particles.set_particles(pos,vel,acc)

        return particles

    def _advance_particles_RK2(self, dt, particles):

        # TODO
        nparticles = particles.nparticles
        mass = particles.masses

        pos = particles.positions
        vel = particles.velocities
        acc = self._calculate_acceleration(nparticles, mass, pos)
    

        # do the RK2 update
        pos2 = pos + vel*dt 
        vel2 = vel + acc*dt
        acc2 = self._calculate_acceleration(nparticles, mass, pos2) 

        pos2 = pos2 + vel2*dt
        vel2 = vel2 + acc2*dt

        # average
        pos = 0.5*(pos + pos2)
        vel = 0.5*(vel + vel2)
        acc = self._calculate_acceleration(nparticles, mass, pos)

        # update the particles
        particles.set_particles(pos, vel, acc)

        
        return particles

    def _advance_particles_RK4(self, dt, particles):
        
        #TODO








        return particles



if __name__ == "__main__":
    
    pass