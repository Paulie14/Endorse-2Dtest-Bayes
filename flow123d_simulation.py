import os
import subprocess
import numpy as np
import itertools
import collections
import shutil
import csv
import ruamel.yaml as yaml
from typing import List

# from bgem.gmsh import gmsh
# from bgem.gmsh import gmsh_io
# from bgem.gmsh import options as gmsh_options
# from bgem.gmsh import heal_mesh

import matplotlib.pyplot as plt
import aux_functions


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


class endorse_2Dtest():

    def __init__(self, config, clean):

        # TODO: set work dir
        self.work_dir = config["work_dir"]
        self.clean = clean
        self._config = config
        self.sample_dir = ""
        self.sample_counter = -1

    def set_parameters(self, data_par):
        param_list = self._config["mcmc_parameters"]
        assert(len(data_par) == len(param_list))

        for idx, param in enumerate(param_list):
            assert(param in self._config["hm_params"])
            self._config["hm_params"][param] = data_par[idx]

    def get_observations(self):
        try:
            res = self.calculate(self._config)
            return res
        except ValueError:
            return [-1000, []]

    def calculate(self, config_dict):
        """
        The program changes to <work_dir> directory.
        does all the data preparation, passing
        running simulation
        extracting results
        """

        # create sample dir
        self.sample_counter = self.sample_counter + 1
        self.sample_dir = os.path.join(config_dict["work_dir"],
                                       "solver_" + str(config_dict["solver_id"]).zfill(2) +
                                       "_sample_" + str(self.sample_counter).zfill(3))
        os.makedirs(self.sample_dir, mode=0o775, exist_ok=True)
        os.chdir(self.sample_dir)

        print("=========================== RUNNING CALCULATION " +
              "solver {} ".format(config_dict["solver_id"]).zfill(2) +
              "sample {} ===========================".format(self.sample_counter).zfill(3))
        print(self.sample_dir)

        # collect only
        if config_dict["collect_only"]:
            return 2, self.collect_results(config_dict)

        print("Creating mesh...")
        # comp_mesh = self.prepare_mesh(config_dict, cut_tunnel=False)
        comp_mesh = self.prepare_mesh(config_dict, cut_tunnel=True)

        mesh_bn = os.path.basename(comp_mesh)
        config_dict["hm_params"]["mesh"] = mesh_bn

        # endorse_2Dtest.read_physical_names(config_dict, comp_mesh)
        print("Creating mesh...finished")

        if config_dict["mesh_only"]:
            return -10, []  # tag, value_list

        # endorse_2Dtest.prepare_hm_input(config_dict)
        print("Running Flow123d - HM...")
        hm_succeed = self.call_flow(config_dict, 'hm_params', result_files=["flow_observe.yaml"])
        if not hm_succeed:
            # raise Exception("HM model failed.")
            # "Flow123d failed (wrong input or solver diverged)"
            return -1, []  # tag, value_list
        print("Running Flow123d - HM...finished")

        if self._config["make_plots"]:
            self.observe_time_plot(config_dict)

        print("Finished computation")

        collected_values = self.collect_results(config_dict)
        return 1, collected_values  # tag, value_list

    # def check_data(self, data, minimum, maximum):
    #     n_times = len(endorse_2Dtest.result_format()[0].times)
    #     if len(data) != n_times:
    #         raise Exception("Data not corresponding with time axis.")
    #
    #     if np.isnan(np.sum(data)):
    #         raise Exception("NaN present in extracted data.")
    #
    #     min = np.amin(data)
    #     if min < minimum:
    #         raise Exception("Data out of given range [min].")
    #     max = np.amax(data)
    #     if max > maximum:
    #         raise Exception("Data out of given range [max].")

    def collect_results(self, config_dict):
        output_dir = config_dict["hm_params"]["output_dir"]

        # the times defined in input
        times = np.array(generate_time_axis(config_dict))
        with open(os.path.join(output_dir, "flow_observe.yaml"), "r") as f:
            loaded_yaml = yaml.load(f, yaml.CSafeLoader)
            points = loaded_yaml['points']
            point_names = [p["name"] for p in points]
            print("Collecting results for observe points: ", point_names)
            data = loaded_yaml['data']
            values = np.array([d["pressure_p0"] for d in data])
            obs_times = np.array([d["time"] for d in data]).transpose()

            # check that observe data are computed at all times of defined time axis
            all_times_computed = np.alltrue(np.isin(times, obs_times))
            if not all_times_computed:
                raise Exception("Observe data not computed at all times as defined by input!")
            # skip the times not specified in input
            t_indices = np.isin(obs_times, times).nonzero()
            values = values[t_indices].transpose()

        if config_dict["clean_sample_dir"]:
            shutil.rmtree(self.sample_dir)

        # flatten to format: [Point0_all_all_times, Point1_all_all_times, Point2_all_all_times, ...]
        res = values.flatten()
        return res

    def call_flow(self, config_dict, param_key, result_files):
        """
        Redirect sstdout and sterr, return true on succesfull run.
        :param result_files: Files to be computed - skip computation if already exist.
        :param param_key: config dict parameters key
        :param config_dict:
        :return:
        """

        params = config_dict[param_key]
        fname = params["in_file"]
        arguments = config_dict["_aux_flow_path"].copy()
        output_dir = "output_" + fname
        config_dict[param_key]["output_dir"] = output_dir
        if all([os.path.isfile(os.path.join(output_dir, f)) for f in result_files]):
            status = True
        else:
            aux_functions.substitute_placeholders(
                os.path.join(config_dict["common_files_dir"],
                             fname + '_tmpl.yaml'),
                fname + '.yaml',
                params)
            # arguments.extend(['--output_dir', output_dir, fname + ".yaml"])
            if not config_dict["run_on_metacentrum"]:
                arguments.extend(['--output_dir', output_dir, fname + ".yaml"])
                print("Running: ", " ".join(arguments))
                with open(fname + "_stdout", "w") as stdout:
                    with open(fname + "_stderr", "w") as stderr:
                        completed = subprocess.run(arguments, stdout=stdout, stderr=stderr)
                print("Exit status: ", completed.returncode)
                status = completed.returncode == 0
            else:
                from mpi4py import MPI

                # arguments.extend(['--output_dir', os.path.abspath(output_dir), os.path.abspath(fname + ".yaml")])
                # print("Running: ", " ".join(arguments))
                # sub_comm = MPI.COMM_SELF.Spawn('flow123d', args=arguments[1:], maxprocs=1)

                # wrap inside bash script to simulate redirection of stderr and stdout from spawned process
                arguments.extend(['--output_dir', output_dir, fname + ".yaml"])
                lines = [
                    '#!/bin/bash',
                    'cd ' + self.sample_dir,
                    " ".join(arguments) + " 2> stderr.log" + " 1> stdout.log"
                ]
                run_script = "flow123d.sh"
                run_script = os.path.abspath(run_script)
                with open(run_script, 'w') as f:
                    f.write('\n'.join(lines))
                sub_comm = MPI.COMM_SELF.Spawn('bash', args=[run_script], maxprocs=1)
                # sub_comm.Barrier()
                # sub_comm.Disconnect()
                status = True

        conv_check = aux_functions.check_conv_reasons(os.path.join(output_dir, "flow123.0.log"))
        print("converged: ", conv_check)
        return status and (conv_check >= 0)







    def prepare_mesh(self, config_dict, cut_tunnel):
        mesh_name = config_dict["mesh_name"]
        if cut_tunnel:
            mesh_name = mesh_name + "_cut"
        mesh_file = mesh_name + ".msh"
        mesh_healed = mesh_name + "_healed.msh"

        if os.path.isfile(os.path.join(config_dict["common_files_dir"], mesh_healed)):
            shutil.copyfile(os.path.join(config_dict["common_files_dir"], mesh_healed), mesh_healed)
            return mesh_healed

        # if not os.path.isfile(mesh_healed):
        #     if not os.path.isfile(mesh_file):
        #         self.make_mesh_A(config_dict, mesh_name, mesh_file, cut_tunnel=cut_tunnel)
        #         # self.make_mesh_B(config_dict, mesh_name, mesh_file, cut_tunnel=cut_tunnel)
        #     hm = heal_mesh.HealMesh.read_mesh(mesh_file, node_tol=1e-4)
        #     hm.heal_mesh(gamma_tol=0.01)
        #     hm.stats_to_yaml(mesh_name + "_heal_stats.yaml")
        #     hm.write()
        #     assert hm.healed_mesh_name == mesh_healed
        return mesh_healed

    # def make_mesh_A(self, config_dict, mesh_name, mesh_file, cut_tunnel):
    #     geom = config_dict["geometry"]
    #     tunnel_mesh_step = geom['tunnel_mesh_step']
    #     dimensions = geom["box_dimensions"]
    #     tunnel_dims = np.array([geom["tunnel_dimX"], geom["tunnel_dimY"]])/2
    #     tunnel_center = geom["tunnel_center"]

    #     print("load gmsh api")
    #     factory = gmsh.GeometryOCC(mesh_name, verbose=True)
    #     gmsh_logger = factory.get_logger()
    #     gmsh_logger.start()
    #     gopt = gmsh_options.Geometry()
    #     gopt.Tolerance = 0.0001
    #     gopt.ToleranceBoolean = 0.001
    #     # gopt.MatchMeshTolerance = 1e-1
    #     gopt.OCCFixSmallEdges = True
    #     gopt.OCCFixSmallFaces = True

    #     # Main box
    #     box = factory.rectangle(dimensions).set_region("box")
    #     side = factory.line([-dimensions[0]/2, 0, 0], [dimensions[0]/2, 0, 0])
    #     sides = dict(
    #         bottom=side.copy().translate([0, -dimensions[1] / 2, 0]),
    #         top   =side.copy().translate([0, +dimensions[1] / 2, 0]),
    #         left  =side.copy().translate([0, +dimensions[0] / 2, 0]).rotate([0, 0, 1], np.pi / 2),
    #         right =side.copy().translate([0, -dimensions[0] / 2, 0]).rotate([0, 0, 1], np.pi / 2)
    #     )

    #     tunnel_disc = factory.disc(tunnel_center, *tunnel_dims)
    #     tunnel_select = tunnel_disc.copy()

    #     print("cutting and fragmenting...")
    #     box_drilled = box.cut(tunnel_disc)
    #     box_fr, tunnel_fr = factory.fragment(box_drilled, tunnel_disc)

    #     print("marking boundary regions...")
    #     box_all = []

    #     b_box_fr = box_fr.get_boundary()
    #     for name, side_tool in sides.items():
    #         isec = b_box_fr.select_by_intersect(side_tool)
    #         box_all.append(isec.modify_regions("." + name))

    #     if cut_tunnel:
    #         b_tunnel_select = tunnel_select.get_boundary()
    #         b_tunnel = b_box_fr.select_by_intersect(b_tunnel_select)
    #         b_tunnel.modify_regions(".tunnel").mesh_step(tunnel_mesh_step)
    #         box_all.extend([box_fr, b_tunnel])
    #     else:
    #         tunnel = tunnel_fr.select_by_intersect(tunnel_select)
    #         tunnel.set_region("tunnel").mesh_step(tunnel_mesh_step)
    #         box_all.extend([box_fr, tunnel])

    #     mesh_groups = [*box_all]

    #     print("meshing...")

    #     factory.keep_only(*mesh_groups)
    #     # factory.remove_duplicate_entities()
    #     factory.write_brep()

    #     min_el_size = tunnel_mesh_step / 2
    #     max_el_size = np.max(dimensions) / 10

    #     mesh = gmsh_options.Mesh()
    #     # mesh.Algorithm = options.Algorithm2d.MeshAdapt # produce some degenerated 2d elements on fracture boundaries ??
    #     # mesh.Algorithm = options.Algorithm2d.Delaunay
    #     # mesh.Algorithm = options.Algorithm2d.FrontalDelaunay
    #     # mesh.Algorithm3D = options.Algorithm3d.Frontal
    #     # mesh.Algorithm3D = options.Algorithm3d.Delaunay

    #     # mesh.Algorithm = gmsh_options.Algorithm2d.FrontalDelaunay
    #     mesh.Algorithm3D = gmsh_options.Algorithm3d.HXT

    #     mesh.ToleranceInitialDelaunay = 0.01
    #     # mesh.ToleranceEdgeLength = fracture_mesh_step / 5
    #     mesh.CharacteristicLengthFromPoints = True
    #     mesh.CharacteristicLengthFromCurvature = True
    #     mesh.CharacteristicLengthExtendFromBoundary = 2
    #     mesh.CharacteristicLengthMin = min_el_size
    #     mesh.CharacteristicLengthMax = max_el_size
    #     mesh.MinimumCirclePoints = 6
    #     mesh.MinimumCurvePoints = 2

    #     # factory.make_mesh(mesh_groups, dim=2)
    #     factory.make_mesh(mesh_groups)

    #     gmsh_log_msgs = gmsh_logger.get()
    #     gmsh_logger.stop()
    #     check_gmsh_log(gmsh_log_msgs)

    #     factory.write_mesh(format=gmsh.MeshFormat.msh2)
    #     os.rename(mesh_name + ".msh2", mesh_file)
    #     # factory.show()

    # def make_mesh_B(self, config_dict, mesh_name, mesh_file, cut_tunnel):
    #     geom = config_dict["geometry"]
    #     tunnel_mesh_step = geom['tunnel_mesh_step']
    #     dimensions = geom["box_dimensions"]
    #     tunnel_dims = np.array([geom["tunnel_dimX"], geom["tunnel_dimY"]]) / 2
    #     tunnel_center = geom["tunnel_center"]

    #     print("load gmsh api")
    #     factory = gmsh.GeometryOCC(mesh_name, verbose=True)
    #     gmsh_logger = factory.get_logger()
    #     gmsh_logger.start()
    #     gopt = gmsh_options.Geometry()
    #     gopt.Tolerance = 0.0001
    #     gopt.ToleranceBoolean = 0.001
    #     # gopt.MatchMeshTolerance = 1e-1
    #     gopt.OCCFixSmallEdges = True
    #     gopt.OCCFixSmallFaces = True

    #     # Main box
    #     box = factory.rectangle(dimensions).set_region("box")
    #     side = factory.line([-dimensions[0] / 2, 0, 0], [dimensions[0] / 2, 0, 0])
    #     sides = dict(
    #         # bottom=side.copy().translate([0, -dimensions[1] / 2, 0]),
    #         top=side.copy().translate([0, +dimensions[1] / 2, 0]),
    #         # left=side.copy().translate([0, +dimensions[0] / 2, 0]).rotate([0, 0, 1], np.pi / 2),
    #         right=side.copy().translate([0, -dimensions[0] / 2, 0]).rotate([0, 0, 1], np.pi / 2)
    #     )

    #     dir_factor = 0.01
    #     d = dir_factor * np.array(dimensions)
    #     dir_sides = dict(
    #         fix_bottom=factory.line([0, 0, 0], [d[0], 0, 0]).translate([-dimensions[0] / 2, -dimensions[1] / 2, 0]),
    #         bottom=factory.line([0, 0, 0], [dimensions[0] - d[0], 0, 0]).translate([-dimensions[0] / 2 + d[0], -dimensions[1] / 2, 0]),
    #         fix_left=factory.line([0, 0, 0], [0, d[1], 0]).translate([-dimensions[0] / 2, -dimensions[1] / 2, 0]),
    #         left = factory.line([0, 0, 0], [0, dimensions[1] - d[1], 0]).translate([-dimensions[0] / 2, -dimensions[1] / 2 + d[1], 0])
    #     )
    #     dir_sides_select = {key: value.copy() for key, value in dir_sides.items()}

    #     tunnel_disc = factory.disc(tunnel_center, *tunnel_dims)
    #     tunnel_select = tunnel_disc.copy()

    #     print("cutting and fragmenting...")
    #     box_drilled = box.cut(tunnel_disc)

    #     dlist = factory.fragment(box_drilled, tunnel_disc, *[v for v in dir_sides.values()])
    #     box_fr = dlist[0]
    #     tunnel_fr = dlist[1]

    #     print("marking boundary regions...")
    #     box_all = []

    #     b_box_fr = box_fr.get_boundary()
    #     for name, side_tool in sides.items():
    #         isec = b_box_fr.select_by_intersect(side_tool)
    #         box_all.append(isec.modify_regions("." + name))
    #     for name, side_tool in dir_sides_select.items():
    #         isec = b_box_fr.select_by_intersect(side_tool)
    #         box_all.append(isec.modify_regions("." + name))

    #     if cut_tunnel:
    #         b_tunnel_select = tunnel_select.get_boundary()
    #         b_tunnel = b_box_fr.select_by_intersect(b_tunnel_select)
    #         b_tunnel.modify_regions(".tunnel").mesh_step(tunnel_mesh_step)
    #         box_all.extend([box_fr, b_tunnel])
    #     else:
    #         tunnel = tunnel_fr.select_by_intersect(tunnel_select)
    #         tunnel.set_region("tunnel").mesh_step(tunnel_mesh_step)
    #         box_all.extend([box_fr, tunnel])

    #     mesh_groups = [*box_all]

    #     print("meshing...")

    #     factory.keep_only(*mesh_groups)
    #     # factory.remove_duplicate_entities()
    #     factory.write_brep()

    #     min_el_size = tunnel_mesh_step / 2
    #     max_el_size = np.max(dimensions) / 10

    #     mesh = gmsh_options.Mesh()
    #     # mesh.Algorithm = options.Algorithm2d.MeshAdapt # produce some degenerated 2d elements on fracture boundaries ??
    #     # mesh.Algorithm = options.Algorithm2d.Delaunay
    #     # mesh.Algorithm = options.Algorithm2d.FrontalDelaunay
    #     # mesh.Algorithm3D = options.Algorithm3d.Frontal
    #     # mesh.Algorithm3D = options.Algorithm3d.Delaunay

    #     # mesh.Algorithm = gmsh_options.Algorithm2d.FrontalDelaunay
    #     mesh.Algorithm3D = gmsh_options.Algorithm3d.HXT

    #     mesh.ToleranceInitialDelaunay = 0.01
    #     # mesh.ToleranceEdgeLength = fracture_mesh_step / 5
    #     mesh.CharacteristicLengthFromPoints = True
    #     mesh.CharacteristicLengthFromCurvature = True
    #     mesh.CharacteristicLengthExtendFromBoundary = 2
    #     mesh.CharacteristicLengthMin = min_el_size
    #     mesh.CharacteristicLengthMax = max_el_size
    #     mesh.MinimumCirclePoints = 6
    #     mesh.MinimumCurvePoints = 2

    #     # factory.make_mesh(mesh_groups, dim=2)
    #     factory.make_mesh(mesh_groups)

    #     gmsh_log_msgs = gmsh_logger.get()
    #     gmsh_logger.stop()
    #     check_gmsh_log(gmsh_log_msgs)

    #     factory.write_mesh(format=gmsh.MeshFormat.msh2)
    #     os.rename(mesh_name + ".msh2", mesh_file)
    #     # factory.show()

    # @staticmethod
    # def prepare_hm_input(config_dict):
    #     """
    #     Prepare FieldFE input file for the TH simulation.
    #     :param config_dict: Parsed config.yaml. see key comments there.
    #     """
    #     bc_pressure_csv = 'bc_pressure_tunnel.csv'
    #     if os.path.exists(bc_pressure_csv):
    #         return
    #
    #     end_time = 17
    #     time_step = 1
    #     times = np.arange(0, end_time, time_step)
    #     n_steps = len(times)
    #     times = np.append(times, end_time)
    #
    #     start_val = 300
    #     end_val = 0
    #     val_step = (end_val-start_val)/n_steps
    #     values = np.arange(start_val, end_val, val_step)
    #     values = np.append(values, end_val)
    #
    #     header = "time bc_pressure_tunnel"
    #     fmt = "%g"
    #     list_rows = np.column_stack((times, values))
    #     np.savetxt(bc_pressure_csv, list_rows, fmt=fmt, delimiter=' ', header=header)
    #     # with open(bc_pressure_csv, 'w', newline='') as csv_file:
    #     #     writer = csv.writer(csv_file)
    #     #     writer.writerow(["time", "bc_pressure_tunnel"])
    #     #     for t, v in zip(times, values):
    #     #         writer.writerow([t, v])

        # PREPARE normal traction on tunnel boundary evolving in time






    # @staticmethod
    # def extract_time_series(yaml_stream, regions, extract):
    #     """
    #
    #     :param yaml_stream:
    #     :param regions:
    #     :return: times list, list: for every region the array of value series
    #     """
    #     data = yaml.load(yaml_stream, yaml.CSafeLoader)['data']
    #     times = set()
    #     reg_series = {reg: [] for reg in regions}
    #
    #     for time_data in data:
    #         region = time_data['region']
    #         if region in reg_series:
    #             times.add(time_data['time'])
    #             power_in_time = extract(time_data)
    #             reg_series[region].append(power_in_time)
    #     times = list(times)
    #     times.sort()
    #     series = [np.array(region_series, dtype=float) for region_series in reg_series.values()]
    #     return np.array(times), np.array(series)
    #
    # @staticmethod
    # def extract_th_results(output_dir, out_regions, bc_regions):
    #     with open(os.path.join(output_dir, "energy_balance.yaml"), "r") as f:
    #         power_times, reg_powers = endorse_2Dtest.extract_time_series(f, bc_regions, extract=lambda frame: frame['data'][0])
    #         power_series = -sum(reg_powers)
    #
    #     with open(os.path.join(output_dir, "Heat_AdvectionDiffusion_region_stat.yaml"), "r") as f:
    #         temp_times, reg_temps = endorse_2Dtest.extract_time_series(f, out_regions, extract=lambda frame: frame['average'][0])
    #     with open(os.path.join(output_dir, "water_balance.yaml"), "r") as f:
    #         flux_times, reg_fluxes = endorse_2Dtest.extract_time_series(f, out_regions, extract=lambda frame: frame['data'][0])
    #     sum_flux = sum(reg_fluxes)
    #
    #     reg_temps = reg_temps - endorse_2Dtest.zero_temperature_offset
    #
    #     avg_temp_flux = sum([temp * flux for temp, flux in zip(reg_temps, reg_fluxes)]) / sum_flux
    #     return avg_temp_flux, power_series

    # @staticmethod
    # def plot_exchanger_evolution(temp_times, avg_temp, power_times, power_series):
    #     year_sec = 60 * 60 * 24 * 365
    #
    #     import matplotlib.pyplot as plt
    #     fig, ax1 = plt.subplots()
    #     temp_color = 'red'
    #     ax1.set_xlabel('time [y]')
    #     ax1.set_ylabel('Temperature [C deg]', color=temp_color)
    #     ax1.plot(temp_times[1:] / year_sec, avg_temp[1:], color=temp_color)
    #     ax1.tick_params(axis='y', labelcolor=temp_color)
    #
    #     ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    #     pow_color = 'blue'
    #     ax2.set_ylabel('Power [MW]', color=pow_color)  # we already handled the x-label with ax1
    #     ax2.plot(power_times[1:] / year_sec, power_series[1:] / 1e6, color=pow_color)
    #     ax2.tick_params(axis='y', labelcolor=pow_color)
    #
    #     fig.tight_layout()  # otherwise the right y-label is slightly clipped
    #     plt.show()

    def observe_time_plot(self, config_dict):

        output_dir = config_dict["hm_params"]["output_dir"]

        with open(os.path.join(output_dir, "flow_observe.yaml"), "r") as f:
            loaded_yaml = yaml.load(f, yaml.CSafeLoader)
            points = loaded_yaml['points']
            point_names = [p["name"] for p in points]
            data = loaded_yaml['data']
            values = np.array([d["pressure_p0"] for d in data]).transpose()
            times = np.array([d["time"] for d in data]).transpose()

            fig, ax1 = plt.subplots()
            temp_color = ['red', 'green', 'violet', 'blue']
            ax1.set_xlabel('time [d]')
            ax1.set_ylabel('pressure [m]')
            for i in range(0, len(point_names)):
                ax1.plot(times, values[i, 0:], color=temp_color[i], label=point_names[i])

            ax1.tick_params(axis='y')
            ax1.legend()

            fig.tight_layout()  # otherwise the right y-label is slightly clipped
            # plt.show()
            plt.savefig("observe_pressure.pdf")
