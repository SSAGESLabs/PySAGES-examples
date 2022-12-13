#!/usr/bin/env python3

"""
CFF simulation of Alanine Dipeptide in vacuum with OpenMM and PySAGES.
"""


# %%
import argparse
import sys
import time

import numpy
import pysages

from pysages.colvars import DihedralAngle
from pysages.methods import CFF
from pysages.utils import try_import

import matplotlib.pyplot as plt
import pickle

openmm = try_import("openmm", "simtk.openmm")
unit = try_import("openmm.unit", "simtk.unit")
app = try_import("openmm.app", "simtk.openmm.app")


# %%
pi = numpy.pi
kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
kB = kB.value_in_unit(unit.kilojoules_per_mole / unit.kelvin)

T = 298.15 * unit.kelvin
dt = 2.0 * unit.femtoseconds
adp_pdb = "adp-vacuum.pdb"

kT = kB*T               # cannot be passed to a `pysages.method` with the units 
kT = 2.4789570296023884 # kJ/mol

# %%
def generate_simulation(pdb_filename=adp_pdb, T=T, dt=dt):
    pdb = app.PDBFile(pdb_filename)

    ff = app.ForceField("amber99sb.xml")
    cutoff_distance = 1.0 * unit.nanometer
    topology = pdb.topology

    system = ff.createSystem(
        topology, constraints=app.HBonds, nonbondedMethod=app.PME, nonbondedCutoff=cutoff_distance
    )

    # Set dispersion correction use.
    forces = {}
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        forces[force.__class__.__name__] = force

    forces["NonbondedForce"].setUseDispersionCorrection(True)
    forces["NonbondedForce"].setEwaldErrorTolerance(1.0e-5)

    positions = pdb.getPositions(asNumpy=True)

    integrator = openmm.LangevinIntegrator(T, 1 / unit.picosecond, dt)
    integrator.setRandomNumberSeed(42)

    # platform = openmm.Platform.getPlatformByName(platform)
    # simulation = app.Simulation(topology, system, integrator, platform)
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy()

    return simulation


# %%
def get_args(argv):
    available_args = [
        ("time-steps", "t", int, 5e5, "Number of simulation steps"),
        ("train-freq", "f", int, 5e3, "Frequency for neural network training"),
    ]
    parser = argparse.ArgumentParser(description="Example script to run CFF")
    for (name, short, T, val, doc) in available_args:
        parser.add_argument("--" + name, "-" + short, type=T, default=T(val), help=doc)
    return parser.parse_args(argv)


# %%
def main(argv=[]):
    args = get_args(argv)

    cvs = [DihedralAngle([4, 6, 8, 14]), DihedralAngle([6, 8, 14, 16])]
    grid = pysages.Grid(lower=(-pi, -pi), upper=(pi, pi), shape=(50, 50), periodic=True)
    topology = (14,)

    # Method
    method = CFF(cvs, grid, topology, kT, train_freq=args.train_freq)

    tic = time.perf_counter()
    raw_result = pysages.run(method, generate_simulation, timesteps = args.time_steps)
    toc = time.perf_counter()
    print(f"Completed the simulation in {toc - tic:0.4f} seconds.")

    # Pickle the results
    pickle.dump( raw_result, open("raw_result.pickle", "wb") )

    # Analysis: Calculate free energy
    result = pysages.analyze(raw_result)
    mesh = result["mesh"]
    A = result["free_energy"]

    # plot and save free energy to a PNG file
    fig, ax = plt.subplots(dpi=120)

    im = ax.imshow(A, interpolation="bicubic", origin="lower", extent=[-pi, pi, -pi, pi])
    ax.contour(A, levels=12, linewidths=0.75, colors="k", extent=[-pi, pi, -pi, pi])
    ax.set_xlabel(r"$\phi$")
    ax.set_ylabel(r"$\psi$")

    cbar = plt.colorbar(im)
    cbar.ax.set_ylabel(r"$A~[kJ/mol]$", rotation=270, labelpad=20)

    fig.savefig("adp-fe.png", dpi=fig.dpi)

    return result


# %%
if __name__ == "__main__":
    main(sys.argv[1:])
