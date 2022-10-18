#!/usr/bin/env python3

"""
SpectralABF simulation of a butane molecule in HOOMD-blue and PySAGES.
"""


# %%
import argparse
import sys
import time

import numpy as np
import hoomd
import hoomd.md
import pysages

from math import sqrt, pi

from pysages import Grid
from pysages.colvars import DihedralAngle
from pysages.methods import SpectralABF

import matplotlib.pyplot as plt


# %%
kT = 0.596161
dt = 0.02045


# %%
def generate_context(kT=kT, dt=dt, **kwargs):
    hoomd.context.initialize("")
    context = hoomd.context.SimulationContext()

    with context:
        snapshot = hoomd.data.make_snapshot(
            N=14,
            box=hoomd.data.boxdim(Lx=41, Ly=41, Lz=41),
            particle_types=["C", "H"],
            bond_types=["CC", "CH"],
            angle_types=["CCC", "CCH", "HCH"],
            dihedral_types=["CCCC", "HCCC", "HCCH"],
            pair_types=["CCCC", "HCCC", "HCCH"],
        )

        snapshot.particles.typeid[0] = 0
        snapshot.particles.typeid[1:4] = 1
        snapshot.particles.typeid[4] = 0
        snapshot.particles.typeid[5:7] = 1
        snapshot.particles.typeid[7] = 0
        snapshot.particles.typeid[8:10] = 1
        snapshot.particles.typeid[10] = 0
        snapshot.particles.typeid[11:14] = 1

        positions = np.array(
            [
                [-2.990196, 0.097881, 0.000091],
                [-2.634894, -0.911406, 0.001002],
                [-2.632173, 0.601251, -0.873601],
                [-4.060195, 0.099327, -0.000736],
                [-2.476854, 0.823942, 1.257436],
                [-2.832157, 1.833228, 1.256526],
                [-2.834877, 0.320572, 2.131128],
                [-0.936856, 0.821861, 1.258628],
                [-0.578833, 1.325231, 0.384935],
                [-0.581553, -0.187426, 1.259538],
                [-0.423514, 1.547922, 2.515972],
                [-0.781537, 1.044552, 3.389664],
                [0.646485, 1.546476, 2.516800],
                [-0.778816, 2.557208, 2.515062],
            ]
        )

        reference_box_low_coords = np.array([-22.206855, -19.677099, -19.241968])
        box_low_coords = np.array(
            [-snapshot.box.Lx / 2, -snapshot.box.Ly / 2, -snapshot.box.Lz / 2]
        )
        positions += box_low_coords - reference_box_low_coords

        snapshot.particles.position[:] = positions[:]

        mC = 12.00
        mH = 1.008

        # fmt: off
        snapshot.particles.mass[:] = [
            mC, mH, mH, mH,  # grouped by carbon atoms
            mC, mH, mH,
            mC, mH, mH,
            mC, mH, mH, mH,
        ]

        reference_charges = np.array(
            [
                -0.180000, 0.060000, 0.060000, 0.060000,  # grouped by carbon atoms
                -0.120000, 0.060000, 0.060000,
                -0.120000, 0.060000, 0.060000,
                -0.180000, 0.060000, 0.060000, 0.060000,
            ]
        )
        # fmt: on

        charge_conversion = 18.22262
        snapshot.particles.charge[:] = charge_conversion * reference_charges[:]

        snapshot.bonds.resize(13)
        snapshot.bonds.typeid[0:3] = 1
        snapshot.bonds.typeid[3] = 0
        snapshot.bonds.typeid[4:6] = 1
        snapshot.bonds.typeid[6] = 0
        snapshot.bonds.typeid[7:9] = 1
        snapshot.bonds.typeid[9] = 0
        snapshot.bonds.typeid[10:13] = 1

        # fmt: off
        snapshot.bonds.group[:] = [
            [0, 2], [0, 1], [0, 3], [0, 4],  # grouped by carbon atoms
            [4, 5], [4, 6], [4, 7],
            [7, 8], [7, 9], [7, 10],
            [10, 11], [10, 12], [10, 13],
        ]
        # fmt: on

        snapshot.angles.resize(24)
        snapshot.angles.typeid[0:2] = 2
        snapshot.angles.typeid[2] = 1
        snapshot.angles.typeid[3] = 2
        snapshot.angles.typeid[4:8] = 1
        snapshot.angles.typeid[8] = 0
        snapshot.angles.typeid[9] = 2
        snapshot.angles.typeid[10:14] = 1
        snapshot.angles.typeid[14] = 0
        snapshot.angles.typeid[15] = 2
        snapshot.angles.typeid[16:21] = 1
        snapshot.angles.typeid[21:24] = 2

        # fmt: off
        snapshot.angles.group[:] = [
            [1, 0, 2], [2, 0, 3], [2, 0, 4],  # grouped by carbon atoms
            [1, 0, 3], [1, 0, 4], [3, 0, 4],
            # ---
            [0, 4, 5], [0, 4, 6], [0, 4, 7],
            [5, 4, 6], [5, 4, 7], [6, 4, 7],
            # ---
            [4, 7, 8], [4, 7, 9], [4, 7, 10],
            [8, 7, 9], [8, 7, 10], [9, 7, 10],
            # ---
            [7, 10, 11], [7, 10, 12], [7, 10, 13],
            [11, 10, 12], [11, 10, 13], [12, 10, 13],
        ]
        # fmt: on

        snapshot.dihedrals.resize(27)
        snapshot.dihedrals.typeid[0:2] = 2
        snapshot.dihedrals.typeid[2] = 1
        snapshot.dihedrals.typeid[3:5] = 2
        snapshot.dihedrals.typeid[5] = 1
        snapshot.dihedrals.typeid[6:8] = 2
        snapshot.dihedrals.typeid[8:11] = 1
        snapshot.dihedrals.typeid[11] = 0
        snapshot.dihedrals.typeid[12:14] = 2
        snapshot.dihedrals.typeid[14] = 1
        snapshot.dihedrals.typeid[15:17] = 2
        snapshot.dihedrals.typeid[17:21] = 1
        snapshot.dihedrals.typeid[21:27] = 2

        # fmt: off
        snapshot.dihedrals.group[:] = [
            [2, 0, 4, 5], [2, 0, 4, 6], [2, 0, 4, 7],  # grouped by pairs of central atoms
            [1, 0, 4, 5], [1, 0, 4, 6], [1, 0, 4, 7],
            [3, 0, 4, 5], [3, 0, 4, 6], [3, 0, 4, 7],
            # ---
            [0, 4, 7, 8], [0, 4, 7, 9], [0, 4, 7, 10],
            [5, 4, 7, 8], [5, 4, 7, 9], [5, 4, 7, 10],
            [6, 4, 7, 8], [6, 4, 7, 9], [6, 4, 7, 10],
            # ---
            [4, 7, 10, 11], [4, 7, 10, 12], [4, 7, 10, 13],
            [8, 7, 10, 11], [8, 7, 10, 12], [8, 7, 10, 13],
            [9, 7, 10, 11], [9, 7, 10, 12], [9, 7, 10, 13],
        ]
        # fmt: on

        snapshot.pairs.resize(27)
        snapshot.pairs.typeid[0:1] = 0
        snapshot.pairs.typeid[1:11] = 1
        snapshot.pairs.typeid[11:27] = 2
        # fmt: off
        snapshot.pairs.group[:] = [
            # CCCC
            [0, 10],
            # HCCC
            [0, 8],
            [0, 9],
            [5, 10], [6, 10],
            [1, 7], [2, 7], [3, 7],
            [11, 4], [12, 4], [13, 4],
            # HCCH
            [1, 5], [1, 6],
            [2, 5], [2, 6],
            [3, 5], [3, 6],
            [5, 8], [6, 8],
            [5, 9], [6, 9],
            [8, 11], [8, 12], [8, 13],
            [9, 11], [9, 12], [9, 13],
        ]
        # fmt: on

        hoomd.init.read_snapshot(snapshot)

        nl_ex = hoomd.md.nlist.cell()
        nl_ex.reset_exclusions(exclusions=["1-2", "1-3", "1-4"])

        lj = hoomd.md.pair.lj(r_cut=12.0, nlist=nl_ex)
        lj.pair_coeff.set("C", "C", epsilon=0.07, sigma=3.55)
        lj.pair_coeff.set("H", "H", epsilon=0.03, sigma=2.42)
        lj.pair_coeff.set("C", "H", epsilon=sqrt(0.07 * 0.03), sigma=sqrt(3.55 * 2.42))

        coulomb = hoomd.md.charge.pppm(hoomd.group.charged(), nlist=nl_ex)
        coulomb.set_params(Nx=64, Ny=64, Nz=64, order=6, rcut=12.0)

        harmonic = hoomd.md.bond.harmonic()
        harmonic.bond_coeff.set("CC", k=2 * 268.0, r0=1.529)
        harmonic.bond_coeff.set("CH", k=2 * 340.0, r0=1.09)

        angle = hoomd.md.angle.harmonic()
        angle.angle_coeff.set("CCC", k=2 * 58.35, t0=112.7 * pi / 180)
        angle.angle_coeff.set("CCH", k=2 * 37.5, t0=110.7 * pi / 180)
        angle.angle_coeff.set("HCH", k=2 * 33.0, t0=107.8 * pi / 180)

        dihedral = hoomd.md.dihedral.opls()
        dihedral.dihedral_coeff.set("CCCC", k1=1.3, k2=-0.05, k3=0.2, k4=0.0)
        dihedral.dihedral_coeff.set("HCCC", k1=0.0, k2=0.0, k3=0.3, k4=0.0)
        dihedral.dihedral_coeff.set("HCCH", k1=0.0, k2=0.0, k3=0.3, k4=0.0)

        lj_special_pairs = hoomd.md.special_pair.lj()
        lj_special_pairs.pair_coeff.set("CCCC", epsilon=0.07, sigma=3.55, r_cut=12.0)
        lj_special_pairs.pair_coeff.set("HCCH", epsilon=0.03, sigma=2.42, r_cut=12.0)
        lj_special_pairs.pair_coeff.set(
            "HCCC", epsilon=sqrt(0.07 * 0.03), sigma=sqrt(3.55 * 2.42), r_cut=12.0
        )

        coulomb_special_pairs = hoomd.md.special_pair.coulomb()
        coulomb_special_pairs.pair_coeff.set("CCCC", alpha=0.5, r_cut=12.0)
        coulomb_special_pairs.pair_coeff.set("HCCC", alpha=0.5, r_cut=12.0)
        coulomb_special_pairs.pair_coeff.set("HCCH", alpha=0.5, r_cut=12.0)

        hoomd.md.integrate.mode_standard(dt=dt)
        integrator = hoomd.md.integrate.nvt(group=hoomd.group.all(), kT=kT, tau=100 * dt)
        integrator.randomize_velocities(seed=42)

    return context


# %%
def get_args(argv):
    available_args = [
        ("use-grids", "g", bool, 1, "Whether to use grid acceleration"),
        ("time-steps", "t", int, 5e5, "Number of simulation steps"),
    ]
    parser = argparse.ArgumentParser(description="Example script to run SpectralABF")
    for (name, short, T, val, doc) in available_args:
        parser.add_argument("--" + name, "-" + short, type=T, default=T(val), help=doc)

    return parser.parse_args(argv)


# %%
def main(argv=[]):
    args = get_args(argv)

    cvs = [DihedralAngle([0, 4, 7, 10])]
    grid = Grid(lower=(-pi,), upper=(pi,), shape=(64,), periodic=True)

    method = SpectralABF(cvs, grid)

    timesteps = args.time_steps
    tic = time.perf_counter()
    run_result = pysages.run(method, generate_context, timesteps)
    toc = time.perf_counter()
    print(f"Completed the simulation in {toc - tic:0.4f} seconds.")

    # extract free energy from result
    result = pysages.analyze(run_result)
    mesh = result["mesh"]
    fes_fn = result["fes_fn"]
    A = fes_fn(mesh)
    A = A.max() - A

    # plot the free energy
    fig, ax = plt.subplots()

    ax.set_xlabel(r"Dihedral Angle, $\xi$")
    ax.set_ylabel(r"$A(\xi)$")

    ax.plot(mesh, A)
    plt.gca()
    fig.savefig("butane-fe.png")

    # write free energy to file
    dih_vs_A = np.stack([mesh[:, 0], A[:, 0]], axis=1)
    np.savetxt("butane-fe.dat", dih_vs_A, header='"dih", "Free energy"')

    return result["free_energy"]


# %%
if __name__ == "__main__":
    main(sys.argv[1:])
