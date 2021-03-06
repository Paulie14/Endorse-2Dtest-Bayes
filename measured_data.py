import os
import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy.fft as fft
from scipy import interpolate
from aux_functions import generate_time_axis

class MeasuredData:
    def __init__(self, config):
        self.measured_data = {}
        self.borehole_names = []
        self.interp_data = {}
        self._config = config

    def initialize(self):
        self.borehole_names, self.measured_data = self.read_chandler_data()

        for bname in self.borehole_names:
            t = self.measured_data[bname]["time"]
            p = self.measured_data[bname]["pressure"]
            # self.interp_data[bname] = interpolate.splrep(t, p)
            self.interp_data[bname] = interpolate.CubicSpline(t, p)

    def generate_measured_samples(self, boreholes):
        times = generate_time_axis(self._config)

        # sample measured data at generated times
        values = []
        for bname in boreholes:
            p = self.interp_data[bname](times)
            # transform pressure from [kPa] to pressure head [m]
            # 1 kPa = [/rho/g] = 0.1 m
            p = p / 10
            values.extend(p)
        return times, values

    def plot_interp_data(self):
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('time [d]')
        ax1.set_ylabel('pressure [m]')

        self.plot_data_set(self.borehole_names, self.measured_data, ax1, linestyle='solid')

        for bname in self.borehole_names:
            t = self.measured_data[bname]["time"]
            # p = interpolate.splev(t, self.interp_data[bname])
            p = self.interp_data[bname](t)
            ax1.plot(t, p, color='black', label=bname, linestyle='dotted')

        ax1.tick_params(axis='y')
        ax1.legend(ncol=3)

        fig.tight_layout()
        # plt.show()

        fig_file = os.path.join(self._config["work_dir"], "interp_data_TSX.pdf")
        plt.savefig(fig_file)

    def plot_all_data(self):
        bnames_ch, data_ch = self.read_chandler_data()
        bnames_rq, data_rq = self.read_rutquist_data()

        fig, ax1 = plt.subplots()
        ax1.set_xlabel('time [d]')
        ax1.set_ylabel('pressure [m]')
        self.plot_data_set(bnames_ch, data_ch, ax1, linestyle='solid')
        self.plot_data_set(bnames_rq, data_rq, ax1, linestyle='dotted')

        ax1.tick_params(axis='y')
        ax1.legend(ncol=3)

        fig.tight_layout()
        # plt.show()
        fig_file = os.path.join(self._config["work_dir"], "measured_data_TSX.pdf")
        plt.savefig(fig_file)

    def read_csv_graph_data(self, csv_file):
        with open(csv_file, newline='') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            line_count = 0
            borehole_names = []
            data = {}
            for row in csv_reader:
                if line_count == 0:
                    for col in row:
                        if col != "":
                            borehole_names.append(col)
                            data[col] = {"time": [], "pressure": []}

                if line_count < 2:
                    line_count = line_count + 1
                    continue
                assert 2*len(borehole_names) == len(row)
                for i in range(0, len(borehole_names)):
                    try:
                        t = float(row[2*i])
                        v = float(row[2*i+1])
                        d = data[borehole_names[i]]
                        d["time"].append(t)
                        d["pressure"].append(v)
                    except:
                        line_count = line_count + 1
                        continue
                line_count = line_count + 1
        return borehole_names, data

    def read_chandler_data(self):
        # read Chandler data
        datafile1 = os.path.join(self._config["script_dir"], "measured_data", "chandler_wpd_datasets_HGT1.csv")
        datafile2 = os.path.join(self._config["script_dir"], "measured_data", "chandler_wpd_datasets_HGT2.csv")
        borehole_names, data = self.read_csv_graph_data(datafile1)
        borehole_names_2, data_2 = self.read_csv_graph_data(datafile2)

        borehole_names.extend(borehole_names_2)
        data.update(data_2)
        excavation_start = 12.7
        # sorting and cropping data
        for bname, dat in data.items():
            t = np.array(dat["time"])
            v = np.array(dat["pressure"])
            permutation = np.argsort(t)
            start_idx = np.searchsorted(t, excavation_start, sorter=permutation)
            t = t - excavation_start
            dat["time"] = t[permutation][start_idx:]
            dat["pressure"] = v[permutation][start_idx:]

        return borehole_names, data

    def read_rutquist_data(self):
        # read Rutquist data
        datafile = os.path.join(self._config["script_dir"], "measured_data", "rutquist_wpd_datasets.csv")
        borehole_names, data= self.read_csv_graph_data(datafile)
        excavation_start_rq = 14
        # sorting and cropping data
        for bname, dat in data.items():
            t = np.array(dat["time"])
            v = np.array(dat["pressure"])
            permutation = np.argsort(t)
            start_idx = np.searchsorted(t, excavation_start_rq, sorter=permutation)
            t = t - excavation_start_rq
            dat["time"] = t[permutation][start_idx:]
            dat["pressure"] = v[permutation][start_idx:]
        return borehole_names, data

    def plot_data_set(self, bnames, data, axes, linestyle):
        temp_color = {'HGT1-1': 'sienna', 'HGT1-2': 'yellow', 'HGT1-3': 'orange', 'HGT1-4': 'green', 'HGT1-5': 'red',
                      'HGT2-1': 'teal', 'HGT2-2': 'cyan', 'HGT2-3': 'blue', 'HGT2-4': 'violet'}
        for bname in bnames:
            axes.plot(data[bname]["time"], data[bname]["pressure"], color=temp_color[bname],
                      label=bname, linestyle=linestyle)


if __name__ == "__main__":

    import json
    import flow_wrapper

    config_dict = flow_wrapper.setup_config()

    os.makedirs(config_dict["work_dir"], mode=0o775, exist_ok=True)
    os.chdir(config_dict["work_dir"])

    md = MeasuredData(config_dict)
    md.initialize()

    md.plot_all_data()
    md.plot_interp_data()

    boreholes = ["HGT1-5", "HGT1-4", "HGT2-4", "HGT2-3"]
    times, values = md.generate_measured_samples(boreholes)

    with open("/home/paulie/Workspace/Endorse-2Dtest-Bayes/minimal_flow_PE.json") as f:
        conf = json.load(f)
        conf["problem_parameters"]["observations"] = values

    with open("/home/paulie/Workspace/Endorse-2Dtest-Bayes/minimal_flow_PE.json", "w") as f:
        json.dump(conf, f, indent=4)



# # Number of samples in normalized_tone
# N = SAMPLE_RATE * DURATION
#
# yf = fft.rfft(normalized_tone)
# xf = fft.rfftfreq(N, 1 / SAMPLE_RATE)
#
# plt.plot(xf, np.abs(yf))
# plt.show()
