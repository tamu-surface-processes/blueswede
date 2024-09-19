import sys
import os

from netCDF4 import Dataset


def open_sww(sww_file):
    """Open an sww file and return a netcdf4 Dataset.

    Parameters
    ----------
    sww_file : str
        File to open.

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
    return swwDict
