### make a grid and gridding function ###
import numpy as np
import matplotlib.pyplot as plt

import scipy.spatial as spatial

import time as time_lib

from datetime import timedelta

import sys
import os

import warnings

import shared

import matplotlib.tri as tri


class InteractiveInspector(object):
    """
    Visualization with interactivity for inspecting the output of an
    ANUGA .sww file.

    Basic user interface:
    * left and right arrows go forward and back in time
    * press "d" to show depth, "h" to show stage, "v" to show velocity, t" to show calculated timestep

    .. important::
        This feature is considered experimental and unsupported. It is
        mostly undocumented. It has very limited flexibiltiy.

    Below is a list of suggested improvements for this to be a true feature.
    - [ ] Colorbar limits need to be determined from the data and update dynamically
    - [ ] Plot aspect ratio should be determined from the domain aspect ratio
    - [ ] Make display options selectable on the screen with radio buttons
    - [ ] Set up as a command line utility with better argument parsing

    Examples
    --------

    You can run this file as a script:

    .. code::
        python ./blueswede/visualization.py --sww_file ./tests/test_data/channel3.sww
    """

    def __init__(self, file):
        """Init.

        Initialize the plot.

        Opens the .sww file at :obj:`file` and shows the first index in the file.

        Parameters
        ----------
        file
            Path to .sww file to open.
        """
        self.FILE = file
        self.IDX = 0

        self.DATA = shared.open_sww(self.FILE)

        X_mesh = self.DATA["x"]
        Y_mesh = self.DATA["y"]
        Z_mesh = self.DATA["elevation"]

        # get the grid data
        x = X_mesh[:].data.astype(float)
        y = Y_mesh[:].data.astype(float)

        t = self.DATA["time"][:].data.astype(float)
        nt = len(t)

        xllcorner = float(self.DATA.xllcorner)
        yllcorner = float(self.DATA.yllcorner)

        # set up figure
        self.FIG, self.AX = plt.subplots(figsize=(6, 6))
        # self.AX = self.AX.flatten()
        plt.subplots_adjust(left=0.085, bottom=0.1, top=0.95, right=0.9)

        # connect keys and mouse
        self.kid = self.FIG.canvas.mpl_connect("key_press_event", self._key_press)

        self.x = np.array(self.DATA.variables["x"])
        self.y = np.array(self.DATA.variables["y"])
        self.triangles = np.array(self.DATA.variables["volumes"])

        vols0 = self.triangles[:, 0]
        vols1 = self.triangles[:, 1]
        vols2 = self.triangles[:, 2]

        self.triang = tri.Triangulation(self.x, self.y, self.triangles)

        self.xc = (self.x[vols0] + self.x[vols1] + self.x[vols2]) / 3.0
        self.yc = (self.y[vols0] + self.y[vols1] + self.y[vols2]) / 3.0

        self.xllcorner = self.DATA.xllcorner
        self.yllcorner = self.DATA.yllcorner
        self.zone = self.DATA.zone

        self.elev = np.array(self.DATA.variables["elevation_c"])
        self.stage = np.array(self.DATA.variables["stage_c"])
        self.xmom = np.array(self.DATA.variables["xmomentum_c"])
        self.ymom = np.array(self.DATA.variables["ymomentum_c"])

        self.minimum_allowed_depth = mad = 0.05

        self.depth = np.zeros_like(self.stage)
        if len(self.elev.shape) == 2:
            self.depth = self.stage - self.elev
        else:
            for i in range(self.depth.shape[0]):
                self.depth[i, :] = self.stage[i, :] - self.elev

        with np.errstate(invalid="ignore", divide="ignore"):
            self.xvel = np.where(self.depth > mad, self.xmom / self.depth, 0.0)
            self.yvel = np.where(self.depth > mad, self.ymom / self.depth, 0.0)

        self.speed = np.sqrt(self.xvel**2 + self.yvel**2)

        self.speed_depth = self.speed * self.depth

        self.time = np.array(self.DATA.variables["time"])
        self.nt = self.time.shape[0]

        # title formatter
        self.new_title = lambda: f"{self.IDX}, {self.var_name}"

        # need to compute size of each element
        x0 = self.x[vols0]
        y0 = self.y[vols0]
        x1 = self.x[vols1]
        y1 = self.y[vols1]
        x2 = self.x[vols2]
        y2 = self.y[vols2]
        self.areas = np.zeros((len(self.triangles)))
        self.radii = np.zeros((len(self.triangles)))
        self.areas[:] = (
            -((x1 * y0 - x0 * y1) + (x2 * y1 - x1 * y2) + (x0 * y2 - x2 * y0)) / 2.0
        )
        for i in range(len(self.triangles)):
            # old_rad = self.radii[i]
            a = np.sqrt((x0[i] - x1[i]) ** 2 + (y0[i] - y1[i]) ** 2)
            b = np.sqrt((x1[i] - x2[i]) ** 2 + (y1[i] - y2[i]) ** 2)
            c = np.sqrt((x2[i] - x0[i]) ** 2 + (y2[i] - y0[i]) ** 2)
            self.radii[i] = self.areas[i] / (2 * (a + b + c))

        with np.errstate(invalid="ignore", divide="ignore"):
            self.deltat = np.where(
                self.depth > mad,
                self.radii / (self.speed),
                0.0,
            )
            self.deltat = np.where(
                self.depth > mad,
                (self.speed) / self.radii,
                0.0,
            )
        self.varset = {
            "depth": {"cmap": "Blues", "vmin": 0, "vmax": 5},
            "stage": {"cmap": "Blues", "vmin": 0, "vmax": 5},
            "speed": {"cmap": "plasma", "vmin": 0, "vmax": 2},
            "deltat": {
                "cmap": "plasma",
                "vmin": 0,
                "vmax": 80,
            },
        }

        # find out if elevation changes with time, if so, need to implement routine to update this field with
        if self.elev.ndim > 1:
            self.elev_fixed = False
            raise NotImplementedError(
                "Routines to update bed with timestepping not implemented."
            )
        else:
            self.elev_fixed = True
            elev_init = self.elev[:]

        self.var_name = "depth"
        depth_init = np.copy(self.depth[self.nt - 1, :])

        # self.triang.set_mask(depth_init > mad)
        self.elev_tri = self.AX.tripcolor(
            self.triang, facecolors=elev_init, cmap="Greys_r"
        )

        # self.triang.set_mask(depth_init < mad)
        new_depth = np.copy(depth_init)
        new_depth[depth_init < mad] = np.nan
        self.var_tri = self.AX.tripcolor(
            self.triang,
            facecolors=new_depth,
            cmap=self.varset[self.var_name]["cmap"],
            vmin=self.varset[self.var_name]["vmin"],
            vmax=self.varset[self.var_name]["vmax"],
        )
        self.FIG.colorbar(self.var_tri, ax=self.AX)

        # for axi in self.AX:
        self.AX.set_xlim(np.min(self.x), np.max(self.x))
        self.AX.set_ylim(np.min(self.y), np.max(self.y))
        self.AX.set_aspect("equal")

        self.AX.set_title(self.new_title())

        # min_deltat = np.min(self.deltat, axis=0)
        # min_deltat[np.all(self.deltat == 0, axis=0)] = np.nan

    def _change_disp_var(self, key):
        if key == "d":
            self.var_name = "depth"
        elif key == "v":
            self.var_name = "speed"
        elif key == "h":
            self.var_name = "stage"
        elif key == "t":
            self.var_name = "deltat"
        else:
            pass  # do nothing

        self.var_tri.set_cmap(self.varset[self.var_name]["cmap"])
        self.var_tri.set_clim(
            self.varset[self.var_name]["vmin"], self.varset[self.var_name]["vmax"]
        )
        self._change_time_idx(
            0
        )  # update the display with a new field without changing time

    def _change_time_idx(self, interval):
        next_idx = self.IDX + interval
        # limit to the range of valid indices
        next_idx = min(next_idx, self.nt - 1)
        self.IDX = max(next_idx, 0)

        # change the variable as needed
        new_depth = self.depth[self.IDX, :]
        if self.var_name == "depth":
            new_data = np.copy(self.depth[self.IDX, :])
        elif self.var_name == "stage":
            new_data = self.stage[self.IDX, :]
        elif self.var_name == "speed":
            new_data = self.speed[self.IDX, :]
        elif self.var_name == "deltat":
            new_data = self.deltat[self.IDX, :]
        else:
            raise RuntimeError

        new_data[new_depth < self.minimum_allowed_depth] = np.nan

        # change the deltat always
        new_deltat = self.deltat[self.IDX, :]
        new_deltat[new_depth < self.minimum_allowed_depth] = np.nan

        self.var_tri.set_array(new_data)
        self.AX.set_title(self.new_title())

        self.FIG.canvas.draw_idle()

    def _key_press(self, event):
        # if self.levee_pick_cnt == 0 and not self.in_levee_pick:
        if event.key == " ":
            self._change_time_idx(interval=1)
        elif event.key == "left":
            self._change_time_idx(interval=-1)
        elif event.key == "right":
            self._change_time_idx(interval=1)
        else:
            # assume it is a variable change
            self._change_disp_var(event.key)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--sww_file", help="input .sww file", type=str)

    args = parser.parse_args()

    if args.sww_file is None:
        raise ValueError("Must specify sww file.")

    app = InteractiveInspector(args.sww_file)
    plt.show()
