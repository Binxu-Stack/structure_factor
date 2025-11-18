#!/usr/bin/env python


import numpy as np

import matplotlib.pyplot as plt
from ase.io import read, write
from dynasor.qpoints import get_spherical_qpoints
from dynasor import compute_static_structure_factors, Trajectory
from dynasor.post_processing import get_spherically_averaged_sample_smearing
import click

@click.command(context_settings=dict(show_default=True))
@click.argument("structure_file", type=str, required=True, default="test/tmp_traj.xyz")
@click.option("-q","--q_max", help="Define the q max.", type=float, default=4.0)
@click.option("-m","--max_points", help="Define the max points.", type=int, default=20000)
@click.option("-w","--q_width", help="Define the q width.", type=float, default=0.02)
@click.option("-n","--nq", help="Define the nq.", type=int, default=1000)
@click.option("-sf","--structure_format", help="Define the format of the structure file.", type=str, default="extxyz")
@click.option("-o","--output_file", help="Define the output file.", type=str, default="sq.dat")
def sq_of_structure_to_file(structure_file:str,
                            q_max:float=4.0,
                            max_points:int=20000,
                            q_width:float=0.02,
                            nq:int=1000,
                            structure_format:str='extxyz',
                            output_file:str=None):
    """
    Calculate the static structure factor S(q) from a structure file.
    Args:
        structure_file: The path to the structure file.
        q_max: The maximum q value.
        max_points: The maximum number of q points.
        q_width: The width of the q points.
        nq: The number of q points.
        output_file: The path to the output file.
        structure_format: The format of the structure file pass to ase.io.read.
    Returns:
        q: The q values.
        Sq: The S(q) values.
    """
    q, Sq = sq_of_structure(structure_file,q_max,max_points,q_width,nq,structure_format)
    with open(output_file,'w') as f:
        for iq, iSq in zip(q, Sq):
            f.write(f"{iq} {iSq}\n")

def sq_of_structure(structure_file:str,
                    q_max:float=4.0,
                    max_points:int=20000,
                    q_width:float=0.02,
                    nq:int=1000,
                    structure_format:str='extxyz'):
    """
    Calculate the static structure factor S(q) from a structure file.
    Args:
        structure_file: The path to the structure file.
        q_max: The maximum q value.
        max_points: The maximum number of q points.
        q_width: The width of the q points.
        nq: The number of q points.
        structure_format: The format of the structure file pass to ase.io.read.
    Returns:
        q: The q values.
        Sq: The S(q) values.
    """
    atoms = read(structure_file, format=structure_format)
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

def calculate_fwhm(x, y):
    """
    Calculate the FWHM of a given x and y values.
    Args:
        x: The x values.
        y: The y values.
    Returns:
        fwhm: The FWHM value.
        xl: The left x value.
        xr: The right x value.
    """
    # 找到峰值
    peak_idx = np.argmax(y)
    peak_value = y[peak_idx]

    # 计算半最大值
    half_max = peak_value / 2
    
    # 找到左右两边半最大值的位置 
    left_idx = np.where(y[:peak_idx] <= half_max)[0][-1]
    right_idx = np.where(y[peak_idx:] <= half_max)[0][0] + peak_idx
    
    # 计算FWHM
    fwhm = x[right_idx] - x[left_idx]
    return fwhm, x[left_idx], x[right_idx]

def fwhm_structure_file(structure_file=None,q_max=4.0,max_points=None,q_width=0.02,nq=1000,structure_format='extxyz',skip_front=100):
    """
    Calculate the FWHM of the static structure factor S(q) from a structure file.
    Args:
        structure_file: The path to the structure file.
        q_max: The maximum q value.
        max_points: The maximum number of q points.
        q_width: The width of the q points for gaussian smearing in q space.
        nq: The number of q points.
        structure_format: The format of the structure file pass to ase.io.read.
        skip_front: The number of points to skip from the front.
    Returns:
        fwhm: The FWHM value in q space in A^-1.
        xl: The left x value in q space in A^-1.
        xr: The right x value in q space in A^-1.
    """
    q, Sq = sq_of_structure(structure_file,q_max=q_max,max_points=max_points,q_width=q_width,nq=nq,structure_format=structure_format)
    fwhm, xl, xr = calculate_fwhm(q[skip_front:],Sq[skip_front:])
    return fwhm, xl, xr

@click.command(context_settings=dict(show_default=True))
@click.argument("structure_file", type=str, required=True, default="test/tmp_traj.xyz")
@click.option("-q","--q_max", help="Define the q max.", type=float, default=4.0)
@click.option("-m","--max_points", help="Define the max points.", type=int, default=40000)
@click.option("-w","--q_width", help="Define the q width.", type=float, default=0.02)
@click.option("-n","--nq", help="Define the nq.", type=int, default=1000)
@click.option("-sf","--structure_format", help="Define the format of the structure file.", type=str, default="extxyz")
@click.option("-s","--skip_front", help="Define the number of points to skip from the front.", type=int, default=100)
@click.option("-o","--output_file", help="Define the output file.", type=str, default="fwhm.dat")
@click.option("-p","--plot", help="Plot the static structure factor S(q) and the FWHM.", is_flag=True, default=False)
@click.option("-sq","--output_sq", help="Whether to output sq.dat", is_flag=True, default=False)
def fwhm_of_structure_to_file(structure_file:str,
                            q_max:float=4.0,
                            max_points:int=20000,
                            q_width:float=0.02,
                            nq:int=1000,
                            structure_format:str='extxyz',
                            skip_front:int=100,
                            output_file:str=None,
                            output_sq:bool=False,
                            plot:bool=False):
    """
    Calculate the FWHM of the static structure factor S(q) from a structure file.
    Args:
        structure_file: The path to the structure file.
        q_max: The maximum q value.
        max_points: The maximum number of q points.
        q_width: The width of the q points in q space.
        nq: The number of q points.
        skip_front: The number of points to skip from the front.
        structure_format: The format of the structure file pass to ase.io.read. e.g. extxyz, lammps-data, lammps-dump-text, cif, etc.
        output_sq: Whether to output sq.dat.
        plot: Whether to plot the static structure factor S(q) and the FWHM.
    Returns:
        fwhm: The FWHM value in q space in A^-1.
        xl: The left x value in q space in A^-1.
        xr: The right x value in q space in A^-1.
    """
    q, Sq = sq_of_structure(structure_file,q_max=q_max,max_points=max_points,q_width=q_width,nq=nq,structure_format=structure_format)
    fwhm, xl, xr = calculate_fwhm(q[skip_front:],Sq[skip_front:])
    print(f"FWHM: {fwhm} A^-1")
    print(f"Left x: {xl} A^-1")
    print(f"Right x: {xr} A^-1")
    with open(output_file,'w') as f:
        f.write(f"FWHM: {fwhm} A^-1\n")
        f.write(f"Left x: {xl} A^-1\n")
        f.write(f"Right x: {xr} A^-1\n")
    if output_sq:
        with open(f"{output_file}.sq.dat",'w') as f:
            for iq, iSq in zip(q[skip_front:], Sq[skip_front:]):
                f.write(f"{iq} {iSq}\n")
    if plot:
        plt.plot(q[skip_front:], Sq[skip_front:])
        plt.axvline(x=xl, color='red')
        plt.axvline(x=xr, color='red')
        plt.xlabel(r"q ($\mathrm{rad}\ \mathrm{\AA}^{-1}$)")
        plt.ylabel("S(q)")
        plt.title("Static Structure Factor S(q)")
        plt.savefig(f"{output_file}.svg")
        #plt.show()
        plt.close()
    return fwhm, xl, xr


if __name__ == "__main__":
    sq_of_structure(structure_file="test/tmp_traj.xyz",structure_format="extxyz")



