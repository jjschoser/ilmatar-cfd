import numpy as np
from pathlib import Path

def get_data_header(fname, lo, hi, res, NVARS, step=0, time=0.0):
    assert len(res) == len(lo) and len(res) == len(hi)
    header = "\n".join([str(step), 
                        str(time), 
                        " ".join([str(x) for x in lo]), 
                        " ".join([str(x) for x in hi]), 
                        " ".join([str(x) for x in res]),
                        str(NVARS),
                        fname + ".dat"])
    return header


def write_data(fname, lo, hi, data, step=0, time=0.0, create_header=True):
    res = np.asarray(data.shape[:-1])
    NVARS = data.shape[-1]
    if create_header:
        header = get_data_header(fname, lo, hi, res, NVARS, step=step, time=time)
        with open(fname + ".txt", "w") as f:
            f.write(header)
    REAL = np.float64  # Assumes that REAL has double precision
    with open(fname + ".dat", "wb") as f:
        buffer = np.ascontiguousarray(data.astype(REAL))
        buffer.tofile(f)


def get_sdf_header(fname, lo, hi, res):
    assert len(res) == len(lo) and len(res) == len(hi)
    header = "\n".join([" ".join([str(x) for x in lo]), 
                        " ".join([str(x) for x in hi]), 
                        " ".join([str(x) for x in res]),
                        fname + ".dat"])
    return header

def write_sdf(fname, lo, hi, sdf, create_header=True):
    if create_header:
        res = np.asarray(sdf.shape) - 2
        header = get_sdf_header(fname, lo, hi, res)
        with open(fname + ".txt", "w") as f:
            f.write(header)
    REAL = np.float64  # Assumes that REAL has double precision
    with open(fname + ".dat", "wb") as f:
        buffer = np.ascontiguousarray(sdf.astype(REAL))
        buffer.tofile(f)


def write_settings(fname, init_fname, final_fname, final_time, lo_bc, hi_bc, gamma=1.4, sdf_fname=None):
    lines = [init_fname + "\n",
             final_fname + "\n",
             str(final_time) + "\n",
             " ".join([str(x) for x in lo_bc]) + "\n",
             " ".join([str(x) for x in hi_bc]) + "\n",
             str(gamma) + "\n"]
    if sdf_fname is not None:
        lines.append(sdf_fname)
    with open(fname + ".txt", "w") as f:
        f.writelines(lines)


def read_data(fname):
    with open(fname + ".txt", "r") as f:
        step = int(f.readline())
        time = float(f.readline())
        lo = [float(i) for i in f.readline().split()]
        hi = [float(i) for i in f.readline().split()]
        res = [int(i) for i in f.readline().split()]
        NVARS = int(f.readline())
        data_fname = f.readline().rstrip()
        assert len(lo) == len(hi) and len(lo) == len(res)
    REAL = np.float64  # Assumes that REAL has double precision
    count = np.prod(res) * NVARS
    with open(fname + ".dat", "rb") as f:
        data = np.fromfile(f, dtype=REAL, count=count).reshape((*res, NVARS))
    return lo, hi, data, step, time


# This function was written with the help of ChatGPT
def get_last_header_filename(name):
    best_num = 0
    best_file = name + "0.txt"

    for path in Path(".").glob(f"{name}*.txt"):
        suffix = path.stem[len(name):]  # part after `name`
        if suffix.isdigit():
            num = int(suffix)
            if num > best_num:
                best_num = num
                best_file = path

    return best_file
