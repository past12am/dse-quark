from os import listdir
from os.path import isfile, join
import typing

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import cm
import matplotlib.ticker as mticker


def sigma_v(p2: np.ndarray, A: np.ndarray, M: np.ndarray):
    return 1.0 / (A * (p2 + np.square(M)))

def sigma_s(p2: np.ndarray, A: np.ndarray, M: np.ndarray):
    return sigma_v(p2, A, M) * M



def build_2D_meshgrid_from_lists(p2_compl_full: np.ndarray, amp_compl_full: np.ndarray, limits: typing.Tuple):
    # TODO currently we assume an order in p2_compl
    sorted_idxs = np.argsort(p2_compl_full)
    p2_compl_full_sorted = p2_compl_full[sorted_idxs]
    amp_compl_full_sorted = amp_compl_full[sorted_idxs]

    real_subgrid_idxs = np.argwhere((limits[0] <= np.real(p2_compl_full_sorted)) & (np.real(p2_compl_full_sorted) <= limits[1])).flatten()

    p2_compl = p2_compl_full_sorted[real_subgrid_idxs]
    amp_compl = amp_compl_full_sorted[real_subgrid_idxs]

    p2_real = np.unique(p2_compl.real)
    p2_imag = np.unique(p2_compl.imag)

    real_size = len(p2_real)
    imag_size = len(p2_imag)

    p2_real_meshgrid = np.zeros((real_size, imag_size))
    p2_imag_meshgrid = np.zeros((real_size, imag_size))
    res_real = np.zeros((real_size, imag_size))
    res_imag = np.zeros((real_size, imag_size))

    real_idx = -1
    imag_idx = -1

    last_p2_imag = -1E20
    last_p2_real = -1E20

    for p2, amp in zip(p2_compl, amp_compl):
        if(np.real(p2) > last_p2_real):
            real_idx += 1
            last_p2_real = np.real(p2)

        if(np.imag(p2) > last_p2_imag):
            imag_idx += 1
            last_p2_imag = np.imag(p2)

        elif(np.imag(p2) < last_p2_imag):
            imag_idx = 0
            last_p2_imag = np.imag(p2)


        p2_real_meshgrid[real_idx, imag_idx] = np.real(p2)
        p2_imag_meshgrid[real_idx, imag_idx] = np.imag(p2)
        res_real[real_idx, imag_idx] = np.real(amp)
        res_imag[real_idx, imag_idx] = np.imag(amp)

    return p2_real, p2_real_meshgrid, p2_imag, p2_imag_meshgrid, res_real, res_imag


def plot_sigma_s_sigma_v_contour(p2_real_meshgrid: np.ndarray,
                                 p2_imag_meshgrid: np.ndarray,
                                 sigma_s_real_meshgrid: np.ndarray,
                                 sigma_s_imag_meshgrid: np.ndarray,
                                 sigma_v_real_meshgrid: np.ndarray,
                                 sigma_v_imag_meshgrid: np.ndarray):
    fig = plt.figure(figsize=(12, 9))

    left = 0.125  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.1   # the bottom of the subplots of the figure
    top = 0.9      # the top of the subplots of the figure
    wspace = 0.25  # the amount of width reserved for blank space between subplots
    hspace = 0.3   # the amount of height reserved for white space between subplots
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)

    # Subplot real sigma_v
    ax = fig.add_subplot(2, 2, 1)
    ax.set_title(f"$\Re(\sigma_v)$")
    ax.contourf(p2_real_meshgrid, p2_imag_meshgrid, sigma_v_real_meshgrid, cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")


    # Subplot imag sigma_v
    ax = fig.add_subplot(2, 2, 2)
    ax.set_title(f"$\Im(\sigma_v)$")
    ax.contourf(p2_real_meshgrid, p2_imag_meshgrid, sigma_v_imag_meshgrid, cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")


    # Subplot real sigma_s
    ax = fig.add_subplot(2, 2, 3)
    ax.set_title(f"$\Re(\sigma_s)$")
    ax.contourf(p2_real_meshgrid, p2_imag_meshgrid, sigma_s_real_meshgrid, cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")


    # Subplot imag sigma_s
    ax = fig.add_subplot(2, 2, 4)
    ax.set_title(f"$\Im(\sigma_s)$")
    ax.contourf(p2_real_meshgrid, p2_imag_meshgrid, sigma_s_imag_meshgrid, cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")

    plt.savefig("./output_complex/plots/sigma_contour.png", dpi=600)

    plt.show()
    plt.close()


def plot_sigma_s_sigma_v_around_zero(p2: np.ndarray, sigma_v: np.ndarray, sigma_s: np.ndarray):
    real_subgrid_idx = np.argwhere(np.abs(np.real(p2)) <= np.abs(np.min(np.real(p2)))).flatten()

    fig = plt.figure(figsize=(10, 9))

    left = 0.125  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.1   # the bottom of the subplots of the figure
    top = 0.9      # the top of the subplots of the figure
    wspace = 0.25  # the amount of width reserved for blank space between subplots
    hspace = 0.3   # the amount of height reserved for white space between subplots
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)

    # Subplot real sigma_v
    ax = fig.add_subplot(2, 2, 1, projection='3d')
    ax.set_title(f"$\Re(\sigma_v)$")
    ax.plot_trisurf(np.real(p2[real_subgrid_idx]), np.imag(p2[real_subgrid_idx]), np.real(sigma_v[real_subgrid_idx]), cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")


    # Subplot imag sigma_v
    ax = fig.add_subplot(2, 2, 2, projection='3d')
    ax.set_title(f"$\Im(\sigma_v)$")
    ax.plot_trisurf(np.real(p2[real_subgrid_idx]), np.imag(p2[real_subgrid_idx]), np.imag(sigma_v[real_subgrid_idx]), cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")


    # Subplot real sigma_s
    ax = fig.add_subplot(2, 2, 3, projection='3d')
    ax.set_title(f"$\Re(\sigma_s)$")
    ax.plot_trisurf(np.real(p2[real_subgrid_idx]), np.imag(p2[real_subgrid_idx]), np.real(sigma_s[real_subgrid_idx]), cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")


    # Subplot imag sigma_s
    ax = fig.add_subplot(2, 2, 4, projection='3d')
    ax.set_title(f"$\Im(\sigma_s)$")
    ax.plot_trisurf(np.real(p2[real_subgrid_idx]), np.imag(p2[real_subgrid_idx]), np.imag(sigma_s[real_subgrid_idx]), cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")

    plt.savefig("./output_complex/plots/sigma_around_zero.png", dpi=600)

    plt.show()
    plt.close()

def plot_sigma_s_sigma_v_positive(p2: np.ndarray, sigma_v: np.ndarray, sigma_s: np.ndarray):
    fig = plt.figure(figsize=(10, 9))

    positive_idxs = np.argwhere(np.real(p2) > 0).flatten()

    # Subplot real sigma_v
    ax = fig.add_subplot(2, 2, 1, projection='3d')
    ax.set_title(f"$\Re(\sigma_v)$")
    ax.plot_trisurf(np.log(np.real(p2[positive_idxs])), np.imag(p2[positive_idxs]), np.real(sigma_v[positive_idxs]), cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda val, pos=None: f"$10^{{{val}}}$"))
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))


    # Subplot imag sigma_v
    ax = fig.add_subplot(2, 2, 2, projection='3d')
    ax.set_title(f"$\Im(\sigma_v)$")
    ax.plot_trisurf(np.log(np.real(p2[positive_idxs])), np.imag(p2[positive_idxs]), np.imag(sigma_v[positive_idxs]), cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda val, pos=None: f"$10^{{{val}}}$"))
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))


    # Subplot real sigma_s
    ax = fig.add_subplot(2, 2, 3, projection='3d')
    ax.set_title(f"$\Re(\sigma_s)$")
    ax.plot_trisurf(np.log(np.real(p2[positive_idxs])), np.imag(p2[positive_idxs]), np.real(sigma_s[positive_idxs]), cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda val, pos=None: f"$10^{{{val}}}$"))
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))


    # Subplot imag sigma_s
    ax = fig.add_subplot(2, 2, 4, projection='3d')
    ax.set_title(f"$\Im(\sigma_s)$")
    ax.plot_trisurf(np.log(np.real(p2[positive_idxs])), np.imag(p2[positive_idxs]), np.imag(sigma_s[positive_idxs]), cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda val, pos=None: f"$10^{{{val}}}$"))
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    plt.savefig("./output_complex/plots/sigma_greater_zero.png", dpi=600)

    plt.show()
    plt.close()

def generatePlots(data_files_M, data_files_A):

    for filename_M, filename_A in zip(data_files_M, data_files_A):
        pd_M_p2 = pd.read_csv(base_path + filename_M, decimal=".", delimiter=";")
        pd_A_p2 = pd.read_csv(base_path + filename_A, decimal=".", delimiter=";")

        mass = filename_M[:-4].split("_")[-1]

        p2 = np.array(pd_M_p2["p2_real"].to_numpy() + 1j * pd_M_p2["p2_imag"], dtype=np.complex128)
        M = np.array(pd_M_p2["M_real"].to_numpy() + 1j * pd_M_p2["M_imag"], dtype=np.complex128)
        A = np.array(pd_A_p2["A_real"].to_numpy() + 1j * pd_A_p2["A_imag"], dtype=np.complex128)

        sigma_s_array = sigma_s(p2=p2, A=A, M=M)
        sigma_v_array = sigma_v(p2=p2, A=A, M=M)

        plot_sigma_s_sigma_v_positive(p2=p2, sigma_s=sigma_s_array, sigma_v=sigma_v_array)
        plot_sigma_s_sigma_v_around_zero(p2=p2, sigma_s=sigma_s_array, sigma_v=sigma_v_array)

        p2_real, p2_real_meshgrid, p2_imag, p2_imag_meshgrid, sigma_s_real_meshgrid, sigma_s_imag_meshgrid = build_2D_meshgrid_from_lists(p2, sigma_s_array, limits=(-2, 2))
        _      , _               , _      , _               , sigma_v_real_meshgrid, sigma_v_imag_meshgrid = build_2D_meshgrid_from_lists(p2, sigma_v_array, limits=(-2, 2))

        plot_sigma_s_sigma_v_contour(p2_real_meshgrid, p2_imag_meshgrid, sigma_s_real_meshgrid, sigma_s_imag_meshgrid, sigma_v_real_meshgrid, sigma_v_imag_meshgrid)




base_path = "./output_complex/"
data_files_m = [f for f in listdir(base_path) if isfile(join(base_path, f)) and f[0] == "M" and f.split(".")[-1] == "txt"]
data_files_a = [f for f in listdir(base_path) if isfile(join(base_path, f)) and f[0] == "A" and f.split(".")[-1] == "txt"]


generatePlots(data_files_m, data_files_a)
