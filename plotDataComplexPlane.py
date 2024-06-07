from os import listdir
from os.path import isfile, join

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import cm


def sigma_v(p2: np.ndarray, A: np.ndarray, M: np.ndarray):
    return 1.0 / (A * (p2 + np.square(M)))

def sigma_s(p2: np.ndarray, A: np.ndarray, M: np.ndarray):
    return sigma_v(p2, A, M) * M


def plot_sigma_s_sigma_v(p2: np.ndarray, sigma_v: np.ndarray, sigma_s: np.ndarray):
    fig = plt.figure(figsize=(10, 9))

    # Subplot real sigma_v
    ax = fig.add_subplot(2, 2, 1, projection='3d')
    ax.set_title(f"$\Re(\sigma_v)$")
    ax.plot_trisurf(np.log(np.real(p2)), np.imag(p2), np.real(sigma_v), cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")


    # Subplot imag sigma_v
    ax = fig.add_subplot(2, 2, 2, projection='3d')
    ax.set_title(f"$\Im(\sigma_v)$")
    ax.plot_trisurf(np.log(np.real(p2)), np.imag(p2), np.imag(sigma_v), cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")


    # Subplot real sigma_s
    ax = fig.add_subplot(2, 2, 3, projection='3d')
    ax.set_title(f"$\Re(\sigma_s)$")
    ax.plot_trisurf(np.log(np.real(p2)), np.imag(p2), np.real(sigma_s), cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")


    # Subplot imag sigma_s
    ax = fig.add_subplot(2, 2, 4, projection='3d')
    ax.set_title(f"$\Im(\sigma_s)$")
    ax.plot_trisurf(np.log(np.real(p2)), np.imag(p2), np.imag(sigma_s), cmap=cm.coolwarm)
    ax.set_xlabel(f"$\Re(p^2)$")
    ax.set_ylabel("$\Im(p^2)$")


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

        plot_sigma_s_sigma_v(p2=p2, sigma_s=sigma_s_array, sigma_v=sigma_v_array)




base_path = "/home/past12am/OuzoCloud/Studium/Physik/6_Semester/SE_Bachelorarbeit/QCD_Intro/QuarkDSE_v3/output_complex/"
data_files_m = [f for f in listdir(base_path) if isfile(join(base_path, f)) and f[0] == "M" and f.split(".")[-1] == "txt"]
data_files_a = [f for f in listdir(base_path) if isfile(join(base_path, f)) and f[0] == "A" and f.split(".")[-1] == "txt"]


generatePlots(data_files_m, data_files_a)
