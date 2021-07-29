import numpy as np


def generate_time_axis(config_dict):
    end_time = float(config_dict["end_time"])
    output_times = config_dict["output_times"]

    # create time axis
    times = []
    for dt in output_times:
        b = float(dt["begin"])
        s = float(dt["step"])
        e = float(dt["end"])
        times.extend(np.arange(b, e, s))
    times.append(end_time)
    return times
