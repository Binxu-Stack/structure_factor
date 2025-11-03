#!/usr/bin/env python


import numpy as np
import sys
import time
import argparse

import matplotlib.pyplot as plt
from ase.io import read, write
from dynasor.qpoints import get_spherical_qpoints
from dynasor import compute_static_structure_factors, Trajectory
from dynasor.post_processing import get_spherically_averaged_sample_smearing
import click

@click.command()
@click.argument("xyz_file", type=str, required=True, default="test/tmp_traj.xyz")
@click.option("-q","--q_max", help="Define the q max.", type=float, default=4.0)
@click.option("-m","--max_points", help="Define the max points.", type=int, default=20000)
@click.option("-w","--q_width", help="Define the q width.", type=float, default=0.02)
@click.option("-n","--nq", help="Define the nq.", type=int, default=1000)
@click.option("-o","--output_file", help="Define the output file.", type=str, default="sq.dat")
def sq_of_xyz_to_file(xyz_file:str,q_max:float=4.0,max_points:int=20000,q_width:float=0.02,nq:int=1000,output_file:str=None):
    q, Sq = sq_of_xyz(xyz_file,q_max,max_points,q_width,nq)
    with open(output_file,'w') as f:
        for iq, iSq in zip(q,Sq):
            f.write(f"{iq} {iSq}\n")



def sq_of_xyz(xyz_file:str,q_max:float=4.0,max_points:int=20000,q_width:float=0.02,nq:int=1000):
    atoms = read(xyz_file)
    atoms.calc = None
    n_atoms = len(atoms)
    traj = [atoms for _ in range(2)]
    write('tmp_traj.xyz',traj)
    traj = Trajectory('tmp_traj.xyz', trajectory_format='extxyz', atomic_indices="read_from_trajectory")

    # q-points
    q_points = get_spherical_qpoints(traj.cell, q_max=q_max,max_points=max_points)
    #q_points = get_spherical_qpoints(traj.cell, q_max=q_max)

    # compute Sq
    sample = compute_static_structure_factors(traj, q_points)
    q_linspace = np.linspace(0,q_max,nq)
    sample_averaged = get_spherically_averaged_sample_smearing(sample, q_norms=q_linspace, q_width=q_width)
    q, Sq = sample_averaged.q_norms, sample_averaged.Sq.reshape(-1)/n_atoms
    return q, Sq


if __name__ == "__main__":
    sq_of_xyz(xyz_file="test/tmp_traj.xyz")



