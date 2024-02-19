#!/usr/bin/env python
import sys
sys.path.append("..")
import numpy as np
import importlib

from droplet import *
import hoomd
import hoomd.md

import pysages
import pickle
import matplotlib.pyplot as plt

from gyana.plot_util import savefig, plot_data_to_str

def plot_energy(energies, centers, epsilons):
    fig, ax = plt.subplots()

    ax.set_xlabel(r"distance R $[\sigma]$")
    ax.set_ylabel(r"free energy $[\epsilon]$")

    # ax.set_ylim((0,850))
    ax.set_xlim((0,6))

    attachement_string = ""
    for i in range(len(energies)):
        center = centers[i]
        free_energy = energies[i]
        free_energy = free_energy
        free_energy -= free_energy[0]
        epsilon = epsilons[i]
        ax.plot(center, free_energy, '.-', label=rf"$\Delta A = {epsilon}$")
        attachement_string += plot_data_to_str((center, free_energy), r"$\Delta\epsilon = {epsilon-1}k_BT$")
    ax.legend(loc="best")
    savefig(fig, attachement_string, "energy.pdf")
    plt.close(fig)

# def plot_hist(result, bins=50):
#     fig, ax = plt.subplots(2, 2)

#     # ax.set_xlabel("CV")
#     # ax.set_ylabel("p(CV)")

#     counter = 0
#     hist_per = len(result["centers"]) // 4 + 1
#     for x in range(2):
#         for y in range(2):
#             for i in range(hist_per):
#                 if counter + i < len(result["centers"]):
#                     center = np.asarray(result["centers"][counter + i])
#                     histo, edges = result["histograms"][counter + i].get_histograms(bins=bins)
#                     edges = np.asarray(edges)[0]
#                     edges = (edges[1:] + edges[:-1]) / 2
#                     ax[x, y].plot(edges, histo, label=f"center {center}")
#                     ax[x, y].legend(loc="best", fontsize="xx-small")
#                     ax[x, y].set_yscale("log")
#             counter += hist_per
#     while counter < len(result["centers"]):
#         center = np.asarray(result["centers"][counter])
#         histo, edges = result["histograms"][counter].get_histograms(bins=bins)
#         edges = np.asarray(edges)[0]
#         edges = (edges[1:] + edges[:-1]) / 2
#         ax[1, 1].plot(edges, histo, label=f"center {center}")
#         counter += 1

#     attachement_string = str(result)
#     savefig(fig, attachement_string, "hist.pdf")
#     plt.close(fig)


def main(argv):
    if len(argv) != 0:
        print("Usage: ./plot_result.py")
        return


    deltaAs = [0.1, 0.2, 0.3, 0.4]
    centers = []
    energies = []
    for delta in deltaAs:
        with open(f"result_{delta}.pkl", "rb") as pkl_file:
            result = pickle.load(pkl_file)
        centers.append(np.asarray(result["centers"]))
        energies.append(np.asarray(result["free_energy"]))

    plot_energy(energies, centers, deltaAs)


if __name__ == "__main__":
    main(sys.argv[1:])
