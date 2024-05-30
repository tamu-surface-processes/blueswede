### make a grid and gridding function ###
import numpy as np
import matplotlib.pyplot as plt

import scipy.spatial as spatial

from netCDF4 import Dataset

import time as time_lib

from datetime import timedelta

import sys
import os

import warnings


def grid_sww_to_netcdf(sww_file, nc_file=None, nc_description=None, dx=10, knn=3):
    """Grid an output (merged) sww file into a netcdf.

    Parameters
    ----------
    sww_file : str
        File to grid.

    nc_file : str, optional
        File to output. If not provided, default is to output a netcdf file in
        the same folder as `sww_file` and with the same basename.

    nc_description : str, optional
        Description entered into NetCDF "description" field. Useful to attach
        some high-level metadata about ANUGA simulation.

    dx : float, optional
        Grid spacing for NetCDF file. Default is 10 m.

    knn : int, optional
        Nunmber of neighbors to use for k-nearest neighbors interpolation.
        Default is 3.

    Examples
    --------

    .. code::

        from blueswede.gridding import grid_sww_to_netcdf
        grid_sww_to_netcdf("/scratch/users", nc_file=None, nc_description=None, dx=10, knn=3)
    """
    filename = os.path.splitext(os.path.basename(sww_file))[0]

    # read ANUGA variables
    swwDict = Dataset(sww_file)
    X_mesh = swwDict["x"]
    Y_mesh = swwDict["y"]
    Z_mesh = swwDict["elevation"]

    # get the grid data
    x = X_mesh[:].data.astype(float)
    y = Y_mesh[:].data.astype(float)

    t = swwDict["time"][:].data.astype(float)
    nt = len(t)

    xllcorner = float(swwDict.xllcorner)
    yllcorner = float(swwDict.yllcorner)

    cellsize = dx  # need to sanitize this to be evenly divisble or anything?

    # Get some dimensions and make x,y grid
    nx = int(np.ceil((max(x) - min(x)) / cellsize) + 1)
    xvect = np.linspace(min(x), min(x) + cellsize * (nx - 1), nx)
    xvect_proj = np.linspace(
        xllcorner + min(x), xllcorner + min(x) + cellsize * (nx - 1), nx
    )  # projected
    ny = int(np.ceil((max(y) - min(y)) / cellsize) + 1)
    yvect = np.linspace(min(y), min(y) + cellsize * (ny - 1), ny)
    yvect_proj = np.linspace(
        yllcorner + min(y), yllcorner + min(y) + cellsize * (ny - 1), ny
    )  # projected
    gridX, gridY = np.meshgrid(xvect, yvect)
    gridX_proj, gridY_proj = np.meshgrid(
        xvect_proj, (yvect_proj)
    )  # coords for netcdf and display

    inputXY = np.array([x[:], y[:]]).transpose()
    gridXY_array = np.array([np.concatenate(gridX), np.concatenate(gridY)]).transpose()
    gridXY_array = np.ascontiguousarray(gridXY_array)

    # Inverse-distance interpolation
    index_qFun = spatial.cKDTree(inputXY)
    NNInfo = index_qFun.query(gridXY_array, k=knn)
    # Weights for interpolation
    nn_wts = 1.0 / (NNInfo[0] + 1.0e-100)
    nn_inds = NNInfo[1]

    def _interp_func(data):
        if isinstance(data, list):
            data = np.array(data)
        denom = 0.0
        num = 0.0
        for i in np.arange(knn):
            denom += nn_wts[:, i]
            num += data[nn_inds[:, i]].astype(float) * nn_wts[:, i]
        gridded_data = num / denom
        gridded_data.shape = (len(yvect), len(xvect))
        return gridded_data

    # initialize the netcdf
    if nc_file is None:
        # if nothing provided, place it in the same place with the same name
        netcdf_path = os.path.join("..", filename + ".nc")
    else:
        # otherwise, set the output nc file according to the input path
        if nc_file[-1] == os.path.sep or os.path.isdir(nc_file):
            # if path is a folder, place it in that folder with the same name
            #   warning: this is not robust to a directory that does not already exist!
            netcdf_path = os.path.join(nc_file, filename + ".nc")
        else:
            # if path is a filename, create it as is
            netcdf_path = nc_file

    print(netcdf_path)
    output_netcdf, _var_list = _initialize_netcdf(
        netcdf_path,
        coordinates=("time", "northing", "easting"),
        dimensions=(nt, ny, nx),
        clobber_netcdf=True,
    )

    # save time, easting, northing
    output_netcdf["easting"][:] = gridX_proj[0, :].astype(float)
    output_netcdf["northing"][:] = gridY_proj[:, 0].astype(float)
    output_netcdf["time"][:] = swwDict["time"][:].astype(float)

    # save elevation (no time component)
    _grid_elev = _interp_func(swwDict["elevation"][:].astype(float))
    output_netcdf["elevation"][:] = _grid_elev

    # save other vars (time, northing, easting)
    save_var_list = list(_var_list.keys())
    save_var_list.remove("elevation")
    for t in np.arange(len(output_netcdf["time"])):
        for _var in save_var_list:
            # grid the field
            _grid = _interp_func(swwDict[_var][t].data.astype(float))

            output_netcdf[_var][t, :, :] = _grid


def _initialize_netcdf(
    netcdf_path, coordinates, dimensions, description=None, clobber_netcdf=True
):
    """
    path :
        Path to output file

    coordinates
        3-tuple of strings with names of coordinates (dim0, dim1, dim2)

    dimennsions
        3-tuple of int with the lengths of each dimension (ndim0, ndim1, ndim2)

    Example
        Most often for gidding ANUGA outputs, this will be time and UTM easting and northing.

        .. code:: initialize_netcdf("path", ("time", "northing", "easting"), (nt, ny, nx))
    """
    # setup a netcdf file with the same name
    _clobber_netcdf = clobber_netcdf

    if (os.path.exists(netcdf_path)) and (_clobber_netcdf is False):
        raise FileExistsError(
            "Existing NetCDF4 output file in target output location: "
            "{file_path}".format(file_path=netcdf_path)
        )
    elif (os.path.exists(netcdf_path)) and (_clobber_netcdf is True):
        _msg = "Replacing existing netCDF file"
        warnings.warn(UserWarning(_msg))
        try:
            output_netcdf.close()
        except:
            pass
        os.remove(netcdf_path)

    output_netcdf = Dataset(netcdf_path, "w", format="NETCDF4")

    if description:
        output_netcdf.description = description
    output_netcdf.history = "Created " + time_lib.ctime(time_lib.time())
    output_netcdf.source = "ANUGA simulation output, gridded"

    # create master dimensions (pulls from `_netcdf_coords`)
    _netcdf_coords = coordinates  # sanitize this?
    ndim0, ndim1, ndim2 = dimensions  # usually nt, ny, nx for ANUGA
    output_netcdf.createDimension(_netcdf_coords[1], ndim1)
    output_netcdf.createDimension(_netcdf_coords[2], ndim2)
    output_netcdf.createDimension(_netcdf_coords[0], ndim0)

    # create master coordinates (as netCDF variables)
    time = output_netcdf.createVariable("time", "f4", (_netcdf_coords[0],))
    time.units = "second"
    # time[:] = swwDict["time"][:].astype(float)

    # new output format is 1d x and y
    easting = output_netcdf.createVariable("easting", "f4", ("easting"))
    northing = output_netcdf.createVariable("northing", "f4", ("northing"))
    # easting[:] = gridX_proj[0, :]
    # northing[:] = gridY_proj[:, 0]

    easting.units = "meter"
    northing.units = "meter"

    # set up variables for output data grids
    def _create_grid_variable(varname, varunits, vartype="f4", vardims=()):
        _v = output_netcdf.createVariable(varname, vartype, vardims)
        _v.units = varunits

    _var_list = dict()
    _var_list["elevation"] = ["elevation", "meters", "f4", _netcdf_coords[1:]]
    _var_list["stage"] = ["stage", "meters", "f4", _netcdf_coords]
    _var_list["xmomentum"] = ["xmomentum", "meters", "f4", _netcdf_coords]
    _var_list["ymomentum"] = ["ymomentum", "meters", "f4", _netcdf_coords]

    for _val in _var_list.keys():
        _create_grid_variable(
            _val,
            _var_list[_val][1],
            _var_list[_val][2],
            _var_list[_val][3],
        )

    return output_netcdf, _var_list


if __name__ == "__main__":
    # test for a single file

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--sww_file", help="input .sww file", type=str)
    parser.add_argument("--nc_file", help="output .nc file", type=str)
    parser.add_argument(
        "--nc_description", help="description written to netcdf file", type=str
    )
    parser.add_argument("--dx", help="gridden cell size", type=float, default=10)
    parser.add_argument("--knn", help="k nearest neighbors", type=int, default=3)
    args = parser.parse_args()

    if args.sww_file is None:
        raise ValueError("Must specify sww file.")

    # sww_file = os.path.join(os.sep, "scratch", "ANUGA", "mcbride_8mmhr_4hr.sww")

    grid_sww_to_netcdf(
        args.sww_file,
        nc_file=args.nc_file,
        nc_description=args.nc_description,
        dx=args.dx,
        knn=args.knn,
    )
