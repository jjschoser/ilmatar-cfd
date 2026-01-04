from file_handler import *


def write_settings(settings_fname, init_header_fname, final_header_fname, final_time, 
                   lo_bc, hi_bc, gamma, out_interval, sdf_header_fname=None):
    lines = [init_header_fname + "\n",
             final_header_fname + "\n",
             str(final_time) + "\n",
             " ".join([str(x) for x in lo_bc]) + "\n",
             " ".join([str(x) for x in hi_bc]) + "\n",
             str(gamma) + "\n",
             str(out_interval) + "\n"]
    if sdf_header_fname is not None:
        lines.append(sdf_header_fname)
    create_dir_for_file(settings_fname)
    with open(settings_fname, "w") as f:
        f.writelines(lines)
