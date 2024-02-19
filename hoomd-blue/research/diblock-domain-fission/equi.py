
import sys
import time
import h5py
import shutil
import numpy as np

import hoomd
hoomd.context.initialize()
import hoomd.hdf5
import hoomd.md


def main(argv):
    with hoomd.context.SimulationContext():
        system = hoomd.init.read_gsd("start.gsd", time_step=0)

        nl = hoomd.md.nlist.tree()
        nl.reset_exclusions([])

        qr = []
        qr += ['temperature','potential_energy','kinetic_energy']
        qr += ["pressure_xx","pressure_yy","pressure_zz","pressure_xy","pressure_xz","pressure_yz"]
        qr += ["bond_harmonic_energy"]
        qr += ["pair_dpd_energy",]

        harmonic = hoomd.md.bond.harmonic()
        harmonic.bond_coeff.set("backbone", k=16/3, r0=0)

        seed = int(time.time())
        dpd = hoomd.md.pair.dpd(r_cut=1.0, nlist=nl, seed=seed, kT=1.)
        dpd.pair_coeff.set("A", "A", A=5., gamma=1.0)
        dpd.pair_coeff.set("A", "B", A=7., gamma=1.0)
        dpd.pair_coeff.set("B", "B", A=5., gamma=1.0)

        hoomd.md.integrate.nve(group=hoomd.group.all())
        hoomd.md.integrate.mode_standard(dt=1e-4)

        # writes out trajectory every 1000 MD steps
        gsd = hoomd.dump.gsd("equi_trajectory.gsd", group=hoomd.group.all(), period=1000, overwrite=True)

        with hoomd.hdf5.File("equi_data.h5", "w") as h5file:
            hoomd.hdf5.log(h5file, quantities=qr, period=1000)
            hoomd.run(1e6, limit_multiple=10000)
        final = hoomd.dump.gsd("equi_final.gsd", group=hoomd.group.all(), period=None, truncate=True, overwrite=True)

        # Copy result of equi to preprocess for run.py
        Nreplica = 14
        for replica_num in range(Nreplica):
            shutil.copy("equi_final.gsd",
                        f"final_{replica_num}.gsd")


if __name__ == "__main__":
    main(sys.argv)
