#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import math


class PlotTestZmpBasedMethodResults(object):
    def __init__(self,
                 method,
                 plot_capture_point=False,
                 plot_comp_time=False,
                 plot_all_time_steps=False,):
        # Setup plot
        self.fig = plt.figure()
        plt.rcParams["font.size"] = 10
        if isinstance(method, list):
            method_list = method
        else:
            method_list = [method]
        axis_list = ["x", "y"]
        subplot_args = [len(method_list), len(axis_list), 0]
        if plot_comp_time:
            subplot_args[0] += 1

        if plot_all_time_steps:
            subplot_args[0] += 1

        # Loop for methods
        self.result_data_list = {}
        for method_idx, method_str in enumerate(method_list):
            # Load file
            result_file_path = "/tmp/Test{}.txt".format(method_str)
            print("[PlotTestZmpBasedMethodResults] Load {}".format(result_file_path))
            self.result_data_list[method_str] = np.genfromtxt(result_file_path, dtype=None, delimiter=None, names=True)
            result_data = self.result_data_list[method_str]

            # Loop for axes
            for axis_idx, axis_str in enumerate(axis_list):
                # Setup subplot
                subplot_args[2] = 2 * method_idx + axis_idx + 1
                ax = self.fig.add_subplot(*subplot_args)

                # Plot
                planned_zmp_key = "planned_zmp_{}".format(axis_str)
                ax.plot(result_data["time"], result_data[planned_zmp_key],
                        color="red", label="planned ZMP")
                ref_zmp_key = "ref_zmp_{}".format(axis_str)
                if ref_zmp_key in result_data.dtype.names:
                    ax.plot(result_data["time"], result_data[ref_zmp_key],
                            color="blue", linestyle="dashed", label="ref ZMP")
                ref_zmp_min_key = "ref_zmp_min_{}".format(axis_str)
                ref_zmp_max_key = "ref_zmp_max_{}".format(axis_str)
                if ref_zmp_min_key in result_data.dtype.names:
                    ax.plot(result_data["time"], result_data[ref_zmp_min_key],
                            color="blue", linestyle="dashed", label="min ZMP")
                    ax.plot(result_data["time"], result_data[ref_zmp_max_key],
                            color="blue", linestyle="dashed", label="max ZMP")
                com_pos_key = "com_pos_{}".format(axis_str)
                ax.plot(result_data["time"], result_data[com_pos_key],
                        color="green", label="CoM")
                capture_point_key = "capture_point_{}".format(axis_str)
                if plot_capture_point and capture_point_key in result_data.dtype.names:
                    ax.plot(result_data["time"], result_data[capture_point_key],
                            color="orange", label="capture point")

                # Set labels, etc.
                ax.set_title("{}-{}".format(method_str, axis_str.upper()))
                ax.set_xlabel("time [s]", labelpad=-2)
                ax.set_ylabel("pos [m]".format(axis_str))
                ax.grid()
                if method_idx == 0 and axis_idx == 0:
                    ax.legend(loc="lower right")

        # Plot computation time
        if plot_comp_time:
            subplot_args[2] = len(method_list) * len(axis_list) + 1
            ax = self.fig.add_subplot(*subplot_args)
            ax.bar(np.arange(len(method_list)),
                   [np.mean(self.result_data_list[method_str]["computation_time"]) for method_str in method_list],
                   yerr=[np.std(self.result_data_list[method_str]["computation_time"]) for method_str in method_list],
                   tick_label=method_list, ecolor="black", capsize=5, align="center", log=True)
            ax.set_title("Computation time")
            plt.xticks(rotation=30, ha="center")
            ax.set_ylabel("time [ms]".format(axis_str))
            ax.grid(axis="y")

        # Plot all time steps
        width_= 0
        step_list = []
        foot_step_list = [[0.0, -0.1],[0.0, 0.1],[0.2, 0.1],[0.4, -0.1],[0.6,0.1],[0.8,-0.1],[0.6,0.1],[0.6,-0.1]]
        select = lambda x: 'Right' if x % 2 == 0 else 'Left'

        #y軸オフセットを計算する
        for i in result_data["time"]:
            if i < 2:
                step_list.append("Both")
            else:
                f, g = math.modf(i)
                if (0 < f and f <= 0.2) or (0.8 < f and f < 1):
                    step_list.append("Both")
                else:
                    step_list.append(select(g))
        
        if plot_all_time_steps:
            # Loop for methods
            self.result_data_list = {}
            for method_idx, method_str in enumerate(method_list):
                # Load file
                result_file_path = "/tmp/Test{}.txt".format(method_str)
                self.result_data_list[method_str] = np.genfromtxt(result_file_path, dtype=None, delimiter=None, names=True)
                result_data = self.result_data_list[method_str]

                # Setup subplot
                subplot_args[2] = len(method_list) * len(axis_list) + plot_comp_time + 1
                ax = self.fig.add_subplot(*subplot_args)

                #Foot step
                for i,j in foot_step_list:
                    if j < 0:
                        r = patches.Rectangle( xy=(i-0.05,j-0.1) , width=0.1, height=0.1)
                    else:
                        r = patches.Rectangle( xy=(i-0.05,j) , width=0.1, height=0.1)
                    ax.add_patch(r)

                ax.plot(result_data["planned_zmp_x"], self.offset(step_list,result_data["planned_zmp_y"]),
                        color="red", label="planned ZMP")
                if "ref_zmp_x" in result_data.dtype.names and "ref_zmp_y" in result_data.dtype.names:
                    ax.plot(result_data["ref_zmp_x"], self.offset(step_list,result_data["ref_zmp_y"]),
                            color="blue", linestyle="dashed", label="ref ZMP")

                # if "ref_zmp_min_x" in result_data.dtype.names and "ref_zmp_min_y" in result_data.dtype.names:
                #     ax.plot(result_data["ref_zmp_min_x"], self.offset(step_list,result_data["ref_zmp_min_y"]),
                #             color="yellow", linestyle="dashed", label="min ZMP")
                    
                # if "ref_zmp_max_x" in result_data.dtype.names and "ref_zmp_max_y" in result_data.dtype.names:
                #     ax.plot(result_data["ref_zmp_max_x"], self.offset(step_list,result_data["ref_zmp_max_y"]),
                #             color="blue", linestyle="dashed", label="max ZMP")
                    
                # ax.plot(result_data["com_pos_x"], self.offset(step_list,result_data["com_pos_y"]),
                #         color="green", label="CoM")
                        
                # if plot_capture_point and "capture_point_x" in result_data.dtype.names and "capture_point_y" in result_data.dtype.names:
                #     ax.plot(result_data["capture_point_x"], self.offset(step_list,result_data["capture_point_y"]),
                #             color="orange", label="capture point")

                # Set labels, etc.
                ax.set_title("{}-{}".format(method_str, "x-y"))
                ax.set_xlabel("pos {}[m]".format("x"))
                ax.set_ylabel("pos {}[m]".format("y"))
                ax.grid()
                if method_idx == 0:
                    ax.legend(loc="lower right")

        # Show
        # plt.tight_layout()
        if len(method_list) == 1:
            self.fig.set_size_inches((12, 4))
            plt.subplots_adjust(left=0.06, bottom=0.14, right=0.98, top=0.90, wspace=0.2, hspace=0.6)
        else:
            self.fig.set_size_inches((12, 13))
            plt.subplots_adjust(left=0.06, bottom=0.10, right=0.98, top=0.96, wspace=0.2, hspace=0.7)
        plt.show()

    def offset(self,step_list,data):
        offset = []
        for i,j in enumerate(data):
            if step_list[i] == "Left":
                offset.append(j+0.1)
            elif step_list[i] == "Right":
                offset.append(j-0.1)
            else:
                offset.append(j)
        return offset


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plot test results of ZMP-based methods.")
    all_methods = [
        "PreviewControlZmp",
        "DdpZmp",
        "DcmTracking",
        "FootGuidedControl",
        "LinearMpcZmp",
        "IntrinsicallyStableMpc",
        "SingularPreviewControlZmp"
    ]
    parser.add_argument("--method", "-m", type=str, default="PreviewControlZmp",
                        choices=all_methods+["All"])
    parser.add_argument("--plot-capture-point", "-pcp", action="store_true")
    parser.add_argument("--plot-comp-time", "-pct", action="store_true")

    args = parser.parse_args()
    if args.method == "All":
        args.method = all_methods

    print(args.plot_capture_point)
    print(args.plot_comp_time)
    plot = PlotTestZmpBasedMethodResults(method=args.method,
                                         plot_capture_point=args.plot_capture_point,
                                         plot_comp_time=args.plot_comp_time,
                                         plot_all_time_steps=True)
