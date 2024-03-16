import numpy
import pysages
from pysages.methods import SpectralABF
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
from pysages.methods import CVRestraints
from pysages.colvars.coordinates import Distance
from ase import units
from ase import Atoms
import matplotlib.pyplot as plt
from ase.io import read, write
from ase.calculators.cp2k import CP2K
import pickle
import os

# System parameters
T = 300
dt = 0.5 * units.fs
friction = 0.05
append = True

class CVLogger:
    def __init__(self, cv_file, log_period):
        self.cv_file = cv_file
        self.log_period = log_period
        self.counter = 0

    def save_cv(self, xi):
        with open(self.cv_file, "a+", encoding="utf8") as f:
            f.write(str(self.counter) + "\t")
            f.write("\t".join(map(str, xi.flatten())) + "\n")

    def __call__(self, snapshot, state, timestep):
        if self.counter >= self.log_period and self.counter % self.log_period == 0:
            self.save_cv(state.xi)

        self.counter += 1

def simulation(T=T, dt=dt, friction=friction):
    cp2k_inp = '''
&FORCE_EVAL
&DFT
    BASIS_SET_FILE_NAME BASIS_PBE
    POTENTIAL_FILE_NAME PBE_POTENTIALS
    &QS
        METHOD GPW
        EXTRAPOLATION PS
        EXTRAPOLATION_ORDER 3
        EPS_DEFAULT 1.0E-10
    &END QS
    &POISSON
        PERIODIC XYZ
    &END POISSON
    &SCF
        EPS_SCF 1.0E-6
        SCF_GUESS ATOMIC
        MAX_SCF 50
        &OT T
            MINIMIZER DIIS
            PRECONDITIONER FULL_SINGLE_INVERSE
        &END OT
    &END SCF
    &XC
        &XC_FUNCTIONAL PBE
        &END XC_FUNCTIONAL
        &VDW_POTENTIAL
            DISPERSION_FUNCTIONAL PAIR_POTENTIAL
            &PAIR_POTENTIAL
            TYPE DFTD3
            PARAMETER_FILE_NAME dftd3.dat
            REFERENCE_FUNCTIONAL PBE
            R_CUTOFF [angstrom] 16.0
            &END PAIR_POTENTIAL
        &END VDW_POTENTIAL
    &END XC
    &LS_SCF
        MAX_SCF 50
    &END LS_SCF
&END DFT
&SUBSYS
    &KIND H
        MASS 2.0
        BASIS_SET DZVP-MOLOPT-PBE-GTH-q1
        POTENTIAL GTH-PBE-q1
    &END KIND
    &KIND O
        BASIS_SET DZVP-MOLOPT-PBE-GTH-q6
        POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND Na
        BASIS_SET DZVP-MOLOPT-PBE-GTH-q9
        POTENTIAL GTH-PBE-q9
    &END KIND
    &KIND Cl
        BASIS_SET DZVP-MOLOPT-PBE-GTH-q7
        POTENTIAL GTH-PBE-q7
    &END KIND
&END SUBSYS
&END FORCE_EVAL
'''

    if os.path.exists('md.traj'):
        traj = Trajectory('md.traj')
        atoms = traj[-1]
    else:
        atoms = read('initial.xyz')
    
    for atom in atoms:
        new_masses = atoms.get_masses()
        if atom.symbol == 'H':
            new_masses[atom.index] = 2.014
            atoms.set_masses(masses=new_masses)

    calc = CP2K(basis_set_file=None,
                potential_file=None,
                basis_set=None,
                pseudo_potential=None,
                xc=None,
                max_scf=None,
                uks=None,
                cutoff=300*units.Ry,
                stress_tensor=False,
                inp=cp2k_inp)

    atoms.set_calculator(calc)
    # MaxwellBoltzmannDistribution(atoms, temperature_K=300)

    dyn = Langevin(atoms,
                    timestep=dt,
                    temperature_K=T,
                    friction=friction,
                    logfile='md.log')

    if os.path.exists('md.traj') and append:
        traj = Trajectory('md.traj', 'a', atoms)
    else:
        traj = Trajectory('md.traj', 'w', atoms)

    dyn.attach(traj.write, interval=1)

    return dyn

# functions for ploting and storing data
def plot_energy(result):
    fig, ax = plt.subplots()

    ax.set_xlabel("CV")
    ax.set_ylabel("Free energy $[\\epsilon]$")

    free_energy = numpy.asarray(result["free_energy"])
    free_energy = free_energy - free_energy.max()
    x = numpy.asarray(result["mesh"])
    ax.plot(x, free_energy, color="teal")

    fig.savefig("energy.png")

def plot_forces(result):
    fig, ax = plt.subplots()

    ax.set_xlabel("CV")
    ax.set_ylabel("Forces $[\\epsilon]$")

    forces = numpy.asarray(result["mean_force"])
    x = numpy.asarray(result["mesh"])
    ax.plot(x, forces, color="teal")

    fig.savefig("forces.png")

def plot_histogram(result):
    fig, ax = plt.subplots()

    ax.set_xlabel("CV")
    ax.set_ylabel("Histogram $[\\epsilon]$")

    hist = numpy.asarray(result["histogram"]) / numpy.nanmax(
        numpy.asarray(result["histogram"])
    )
    x = numpy.asarray(result["mesh"])
    ax.plot(x, hist, color="teal")

    fig.savefig("histogram.png")

def save_energy_forces(result):
    Energy = numpy.asarray(result["free_energy"])
    Forces = numpy.asarray(result["mean_force"])
    Grid = numpy.asarray(result["mesh"])
    hist = numpy.asarray(result["histogram"]) / numpy.nanmax(
        numpy.asarray(result["histogram"])
    )
    numpy.savetxt("FES.csv", numpy.column_stack([Grid, Energy]))
    numpy.savetxt("Forces.csv", numpy.column_stack([Grid, Forces]))
    numpy.savetxt("Histogram.csv", numpy.column_stack([Grid, hist]))

def main():

    # Define CV grid and CV restraints
    grid = pysages.Grid(lower=(2.8), upper=(3.3), shape=(20,), periodic=False)
    restraint = CVRestraints(lower=(2.75), upper=(3.35), ku=(75), kl=(75))
    
    # Define CV and spectral abf method
    cvs = [Distance([[339],[340]]),]
    method = SpectralABF(cvs, grid, restraints=restraint, N=250)

    # Setup CVlogger to print info about CV
    stride = 1
    cv_file = "cv.dat"
    callback = CVLogger(cv_file, stride)
    
    # Run the method and pickle the results
    # with open('1/raw_result.pickle', 'rb') as f:
    #     state = pickle.load(f)

    state = pysages.run(method,simulation,10,callback)
    with open('raw_result.pickle', 'wb') as f:
        pickle.dump(state, f)

    result = pysages.analyze(state)
    plot_energy(result)
    plot_forces(result)
    plot_histogram(result)
    save_energy_forces(result)

main()