import sys
import time
import datetime
import h5py
import numpy as np
import importlib

from droplet import *

import hoomd
import hoomd.md

import pysages
from pysages.colvars import Distance
from pysages.methods import UmbrellaIntegration, SerialExecutor

import gsd
import gsd.hoomd

import pickle

def generate_context(**kwargs):
    # if kwargs.get("mpi", False):
    #     MPI = importlib.import_module("mpi4py.MPI")
    #     init_kwargs = {"mpi_comm": MPI.COMM_SELF}
    # else:
    #     init_kwargs = {}
    # hoomd.context.initialize("--single-mpi", **init_kwargs)
    hoomd.context.initialize("--single-mpi")

    context = hoomd.context.SimulationContext()
    with context:
        deltaA = kwargs.get("deltaA")
        print(f"Operating replica {kwargs.get('replica_num')}: {deltaA}")
        system = hoomd.init.read_gsd(f"final_{kwargs.get('replica_num')}.gsd", time_step=0)
        box = system.box
        nl = hoomd.md.nlist.tree()
        nl.reset_exclusions([])

        # Note that pysages energies from the bias potential can not be logged from hoomd
        qr = []
        qr += ['temperature','potential_energy','kinetic_energy']
        qr += ["bond_harmonic_energy"]
        qr += ["pair_dpd_energy",]

        harmonic = hoomd.md.bond.harmonic()
        harmonic.bond_coeff.set("backbone", k=16/3, r0=0)

        seed = int(time.time()) ^ kwargs.get("replica_num")
        dpd = hoomd.md.pair.dpd(r_cut=1.0, nlist=nl, seed=seed, kT=1.)
        dpd.pair_coeff.set("A", "A", A=5., gamma=1.0)
        dpd.pair_coeff.set("A", "B", A=5.+deltaA, gamma=1.0)
        dpd.pair_coeff.set("B", "B", A=5., gamma=1.0)

        hoomd.md.integrate.nve(group=hoomd.group.all())
        hoomd.md.integrate.mode_standard(dt=1e-3)

        hoomd.analyze.log(f"data_{kwargs.get('replica_num')}_{deltaA}.dat", qr, period=10000, overwrite=True)
    return context

def post_run_action(**kwargs):
    hoomd.util.quiet_status()
    deltaA = kwargs.get("deltaA")
    hoomd.dump.gsd(
        filename=f"final_{kwargs.get('replica_num'):03d}_{deltaA}.gsd",
        overwrite=True,
        period=None,
        group=hoomd.group.all(),
    )
    hoomd.util.unquiet_status()

def get_executor(mpi):
    if mpi:
        futures = importlib.import_module("mpi4py.futures")
        return futures.MPIPoolExecutor()
    return SerialExecutor()


def main(argv):

    # Find indeces of the droplet forming polymer part
    with gsd.hoomd.open("../start.gsd", "rb") as start_file:
        frame = start_file[0]
        idx = np.where(frame.particles.typeid == frame.particles.types.index("A"))[0].astype("int")
        box = frame.configuration.box[:3]

    Nidx = len(idx)
    idxA = list(idx)[:Nidx//2]
    idxB = list(idx)[Nidx//2:]

    start = datetime.datetime.now()

    for deltaA in [0.0, 0.1, 0.2, 0.3, 0.4]:
        cvs = [Distance((idxA, idxB))]
        centers = list(np.linspace(0, 6, 14))
        method = UmbrellaIntegration(cvs, 1000., centers, 100, int(1e6))

        context_args = {"mpi": True}
        context_args["deltaA"] = deltaA

        executor = get_executor(context_args["mpi"])

        raw_result = pysages.run(method, generate_context, 2e6,
                                 context_args=context_args,
                                 post_run_action=post_run_action,
                                 executor=executor)
        result = pysages.analyze(raw_result)
        print("end", time.asctime())
        end = datetime.datetime.now()
        diff = end - start
        with open(f"result_{deltaA}.pkl", "wb") as file_handle:
            pickle.dump(result, file_handle)


if __name__ == "__main__":
    main(sys.argv)
