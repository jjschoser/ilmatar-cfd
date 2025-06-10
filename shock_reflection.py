import argparse
from matplotlib import pyplot as plt
import numpy as np
import subprocess


def write_data(fname, lo, hi, data, step=0, time=0.0):
    res = np.asarray(data.shape[:-1])
    assert len(res) == len(lo) and len(res) == len(hi)
    NVARS = data.shape[-1]
    header = "\n".join([str(step), 
                        str(time), 
                        " ".join([str(x) for x in lo]), 
                        " ".join([str(x) for x in hi]), 
                        " ".join([str(x) for x in res]),
                        str(NVARS)])
    ncells = np.prod(res)
    np.savetxt(fname, data.reshape(ncells, NVARS), header=header, comments="")


def write_sdf(fname, lo, hi, sdf):
    res = np.asarray(sdf.shape) - 2
    assert len(res) == len(lo) and len(res) == len(hi)
    header = "\n".join([" ".join([str(x) for x in lo]), 
                        " ".join([str(x) for x in hi]), 
                        " ".join([str(x) for x in res])])
    np.savetxt(fname, sdf.flatten(), header=header, comments="")


def write_settings(fname, init_data_fname, final_data_fname, final_time, lo_bc, hi_bc, gamma=1.4, sdf_fname=None):
    with open(fname, "w") as f:
        lines = [init_data_fname + "\n",
                 final_data_fname + "\n",
                 str(final_time) + "\n",
                 " ".join([str(x) for x in lo_bc]) + "\n",
                 " ".join([str(x) for x in hi_bc]) + "\n",
                 str(gamma) + "\n"]
        if sdf_fname is not None:
            lines.append(sdf_fname)
        f.writelines(lines)


def read_data(fname):
    with open(fname, "r") as f:
        step = int(f.readline())
        time = float(f.readline())
        lo = [float(i) for i in f.readline().split()]
        hi = [float(i) for i in f.readline().split()]
        res = [int(i) for i in f.readline().split()]
        NVARS = int(f.readline())
        assert len(lo) == len(hi) and len(lo) == len(res)
    return lo, hi, np.loadtxt(fname, skiprows=6).reshape((*res, NVARS)), step, time


if __name__ == "__main__":
    settings_fname = "ShockReflectionSettings.txt"
    init_data_fname = "ShockReflectionInitData.txt"
    final_data_fname = "ShockReflectionFinalData.txt"
    sdf_fname = "ShockReflectionSDF.txt"

    gamma = 1.4
    x_shock = 4e-3
    x_wedge = 4.96e-3
    alpha_wedge = np.deg2rad(25.0)
    final_time = 35e-6
    rho_inf = 1.225
    vel_inf = 0
    p_inf = 101325
    M_shock = 1.7

    lo_bc = [0, 1]
    hi_bc = [0, 0]

    lo = np.array([-4e-3, 0.0])
    hi = np.array([29e-3, 16.5e-3])
    res = np.array([512, 256])
    
    c_inf = np.sqrt(gamma * p_inf / rho_inf)
    s_shock = M_shock * c_inf
    M_inf = vel_inf / c_inf
    rho_star = rho_inf * ((gamma + 1) * (M_inf - M_shock) ** 2) / ((gamma - 1) * (M_inf - M_shock) ** 2 + 2)
    p_star = p_inf * ((2 * gamma * (M_inf - M_shock) ** 2 - (gamma - 1)) / (gamma + 1))
    vel_star = (1 - rho_inf / rho_star) * s_shock + vel_inf * (1 - rho_inf / rho_star)

    NVARS = 4
    init_data = np.zeros((*res, NVARS))
    sdf = np.zeros(res)
    dx = (hi - lo) / res
    x, y = [np.linspace(lo[d] + 0.5 * dx[d], hi[d] - 0.5 * dx[d], res[d]) for d in range(2)]
    X, Y = np.meshgrid(x, y, indexing="ij")

    init_data[:, :, 0] = np.where(X < x_shock, rho_star, rho_inf)
    vel_x = np.where(X < x_shock, vel_star, vel_inf)
    vel_y = np.zeros_like(X)
    p = np.where(X < x_shock, p_star, p_inf)
    init_data[:, :, 1] = init_data[:, :, 0] * vel_x
    init_data[:, :, 2] = init_data[:, :, 0] * vel_y
    init_data[:, :, 3] = 0.5 * init_data[:, :, 0] * (vel_x ** 2 + vel_y ** 2) + p / (gamma - 1)
    write_data(init_data_fname, lo, hi, init_data)

    x_sdf, y_sdf = [np.linspace(lo[d] - 0.5 * dx[d], hi[d] + 0.5 * dx[d], res[d] + 2) for d in range(2)]
    X_sdf, Y_sdf = np.meshgrid(x_sdf, y_sdf, indexing="ij")
    sdf = np.cos(alpha_wedge) * Y_sdf - np.sin(alpha_wedge) * (X_sdf - x_wedge)
    write_sdf(sdf_fname, lo, hi, sdf)

    write_settings(settings_fname, init_data_fname, final_data_fname, final_time, lo_bc, hi_bc, gamma=gamma, sdf_fname=sdf_fname)
    subprocess.run(["make"])
    subprocess.run(["./simple-cfd", settings_fname])

    _, __, final_data, ___, ____ = read_data(final_data_fname)
    rho = np.where(sdf[1:-1, 1:-1] < 0, np.nan, final_data[:, :, 0])
    
    plt.figure(figsize=(10, 5))
    plt.pcolormesh(X, Y, rho)
    plt.colorbar(label="Density")
    plt.contour(X, Y, sdf[1:-1, 1:-1], levels=[0], colors="k")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(final_data_fname.replace(".txt", ".png"), dpi=300)
    plt.close()
