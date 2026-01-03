import numpy as np
from file_handler import *


def save_data(data_fname, data, double=True):
    REAL = np.float64 if double else np.float32
    create_dir_for_file(data_fname)
    with open(data_fname, "wb") as f:
        buffer = np.ascontiguousarray(data.astype(REAL))
        buffer.tofile(f)


def load_data(data_fname, res, NVARS, double=True):
    REAL = np.float64 if double else np.float32
    count = np.prod(res) * NVARS
    with open(data_fname, "rb") as f:
        data = np.fromfile(f, dtype=REAL, count=count).reshape((*res, NVARS))
    return data


def save(header_fname, step, time, lo, hi, data):
    res = np.asarray(data.shape[:-1])
    NVARS = data.shape[-1]
    write_header(header_fname, step, time, lo, hi, res, NVARS)
    save_data(get_data_fname(header_fname), data)


def write_header(header_fname, step, time, lo, hi, res, NVARS):
    assert len(res) == len(lo) and len(res) == len(hi)
    header = "\n".join([str(step), 
                        str(time), 
                        " ".join([str(x) for x in lo]), 
                        " ".join([str(x) for x in hi]), 
                        " ".join([str(x) for x in res]),
                        str(NVARS),
                        remove_path(get_data_fname(header_fname))])
    create_dir_for_file(header_fname)
    with open(header_fname, "w") as f:
            f.write(header)


def load(header_fname):
    with open(header_fname, "r") as f:
        step = int(f.readline())
        time = float(f.readline())
        lo = [float(i) for i in f.readline().split()]
        hi = [float(i) for i in f.readline().split()]
        res = [int(i) for i in f.readline().split()]
        NVARS = int(f.readline())
        data_fname = f.readline().rstrip()
        assert len(lo) == len(hi) and len(lo) == len(res)
    data = load_data(add_path(header_fname, data_fname), res, NVARS)
    return step, time, lo, hi, data


def save_sdf(header_fname, lo, hi, sdf):
    res = np.asarray(sdf.shape) - 2
    write_sdf_header(header_fname, lo, hi, res)
    save_data(get_data_fname(header_fname), sdf)


def write_sdf_header(header_fname, lo, hi, res):
    assert len(res) == len(lo) and len(res) == len(hi)
    header = "\n".join([" ".join([str(x) for x in lo]), 
                        " ".join([str(x) for x in hi]), 
                        " ".join([str(x) for x in res]),
                        remove_path(get_data_fname(header_fname))])
    create_dir_for_file(header_fname)
    with open(header_fname, "w") as f:
            f.write(header)


def load_sdf(header_fname):
    with open(header_fname, "r") as f:
        lo = [float(i) for i in f.readline().split()]
        hi = [float(i) for i in f.readline().split()]
        res = [int(i) for i in f.readline().split()]
        data_fname = f.readline().rstrip()
        assert len(lo) == len(hi) and len(lo) == len(res)
    data = load_data(add_path(header_fname, data_fname), np.asarray(res) + 2, 1)[..., 0]
    return lo, hi, data
