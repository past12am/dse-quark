from os import listdir
from os.path import isfile, join

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plotLogLog(data_files, data_col_y, target_filename_suffix):
    plt.figure(figsize=(8,5))

    for filename in sorted(data_files):
        pd_M_p2 = pd.read_csv(base_path + filename, decimal=".", delimiter=";")

        mass = filename[:-4].split("_")[-1]

        plt.loglog(pd_M_p2["p2"].to_numpy(), pd_M_p2[data_col_y].to_numpy(), label=mass + " GeV")
        plt.xlabel("$p^2$ / GeV$^2$")
        plt.ylabel("$M(p^2)$ / GeV" if data_col_y == "M" else "$A(p^2)$")


    plt.legend()
    plt.savefig("output/plots/" + target_filename_suffix + ".png", dpi=400)
    plt.show()


base_path = "./output/"
data_files_m = [f for f in listdir(base_path) if isfile(join(base_path, f)) and f[0] == "M" and f.split(".")[-1] == "txt"]
data_files_a = [f for f in listdir(base_path) if isfile(join(base_path, f)) and f[0] == "A" and f.split(".")[-1] == "txt"]


plotLogLog(data_files_m, "M", "M")
plotLogLog(data_files_a, "A", "A")
