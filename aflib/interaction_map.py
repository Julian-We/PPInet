from Bio.PDB.PDBParser import PDBParser
import numpy as np
from time import time
import matplotlib.pyplot as plt
from skimage.filters import gaussian as gaussian_filter
import pickle


def shortest_distance_residues(residueA, residueB):
    # Set te minimum distance to infinity
    min_distance = float('inf')

    # Iterate over all the atoms in the residues and calculate the distance between them
    for atomA in residueA:
        for atomB in residueB:
            distance = atomA - atomB
            if distance < min_distance:
                min_distance = distance
    return round(min_distance, 4)  # Round the distance to 4 decimal places


def get_distance_map(pdb_path, plddt_path='', structure_name='arschkipf'):
    # Load the structure from the pdb file
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(structure_name, pdb_path)

    if plddt_path:
        with open(plddt_path, 'rb') as f:
            plddt_data = pickle.load(f)
            plddt_full = plddt_data['plddt']

    # Get the chains from the structure, and raise error is there are more or less than 2 chains
    list_of_chains = list(structure[0].get_chains())
    if len(list_of_chains) == 2:
        chainAid = list_of_chains[0].id
        chainBid = list_of_chains[1].id
    else:
        raise ValueError('The number of chains is not 2')

    # Get the chains from the structure
    chainA = structure[0][chainAid]
    chainB = structure[0][chainBid]

    # Get a list of the residues from the chains
    chainA_residues = list(chainA.get_residues())
    chainB_residues = list(chainB.get_residues())

    if plddt_path:
        plddt_chainA = plddt_full[:len(chainA_residues)]
        plddt_chainB = plddt_full[len(chainA_residues):]
        if not len(plddt_chainA) == len(chainA_residues) or not len(plddt_chainB) == len(chainB_residues):
            raise ValueError('The length of the plddt values is not the same as the length of the residues')

    # Calculate the distance between each residue in chainA and chainB
    distance_array = np.empty((len(chainA_residues), len(chainB_residues)))
    if plddt_path:
        plddt_array = np.empty((len(chainA_residues), len(chainB_residues)))
    for idA, chainA_curr_res in enumerate(chainA_residues):
        for idB, chainB_curr_res in enumerate(chainB_residues):
            distance_array[idA, idB] = shortest_distance_residues(chainA_curr_res, chainB_curr_res)
            if plddt_path:
                plddt_array[idA, idB] = (plddt_chainA[idA] + plddt_chainB[idB]) / 2

    return (distance_array, plddt_array) if plddt_path else distance_array


def plot_distance_map(distance_array,
                      plddt_array=None,
                      axes=None,
                      distance_threshold=8,
                      name_chainA='Chain A',
                      name_chainB='Chain B',
                      binary_map=False,
                      gaussian=False, **kwargs):
    if axes is None:
        fig, axes = plt.subplots()

    if binary_map:
        distance_array_plot = np.zeros(distance_array.shape)
        distance_array_plot[distance_array < distance_threshold] = 1

    else:
        distance_array_plot = distance_array.copy()

    if gaussian:
        distance_array_plot = gaussian_filter(distance_array_plot, **kwargs)
    axes.set_ylabel(name_chainA)
    axes.set_xlabel(name_chainB)
    if not binary_map:
        axes.imshow(distance_array_plot, cmap='viridis_r', vmin=0,
                    vmax=distance_threshold)
    else:
        axes.imshow(distance_array_plot, cmap='viridis', vmin=0,
                    vmax=1)
    if plddt_array is not None:
        if not binary_map:
            print("PLDDT works best with binary_map=True")
        axes.imshow(plddt_array,
                    cmap='magma_r',
                    alpha=(plddt_array / plddt_array.max()*0.5))

