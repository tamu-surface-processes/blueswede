"""
Utilities for preprocessing and preparing spatial data for simulation.
"""
import rasterio
import pandas as pd
import numpy as np
import os

from anuga.geometry.polygon import is_complex


def read_polygon(
    file,
    id_col="shapeid",
    keep_col=["x", "y"],
    drop_looped=True,
    write=False,
    closed=True,
    verbose=False,
    **kwargs
):
    """Process a .csv file with polygon defined to a polygon for ANUGA.

    Originally written to operate on the output from MMQGIS plugin with csv
    created by the "Geometry Export to CSV File" utlitiy, this function has
    been made more flexible by the addition of several input parameters that
    control the processing workflow.

    This function is written as a flexible drop-in replacement for the
    anuga.read_polygon feature.

    Parameters
    ----------
    file : str
        Path to csv file to process.

    keep_col : list
        Which columns to keep for the x-y points, as a list. Default is `
        ["x", "y"]`.

    drop_looped
        Whether to drop a loop in the polygon if one is found.

    write
        Whether to write out the cleaned file to disk.

    closed
        Whether the polygon you are loading is expected to be closed or not.
        This should be `False` for loading polylines, and `True` for loading
        polygons. Default is `True`.

    verbose
        Whether to print certain warnings.

    kwargs
        Passed to pd.read_csv() as kwargs.

    Returns
    -------
    polygons : np.array or list of np.array
        A single numpy array if a single polygon is defined in the csv
        file :obj:`file`. If multiple polygons defined, then a list of
        np.array for each polygon is returned.
    """

    def _parse_and_drop_and_validate(full, uid):
        """private parser.

        Takes full dataframe and unique id.
        """
        # pull out only the points matching the uid
        new_f = f[f[id_col] == uid]
        # keep only the x and y now
        new_f = new_f[keep_col]
        if drop_looped:
            # check first if there is a loop
            looped = np.all(new_f.iloc[0] == new_f.iloc[-1])
            if looped:
                # drop the last item (i.e., the loop)
                new_f = new_f[:-1]

        # make it into a polygon as ANUGA expects (a list of lists)
        polygon = new_f.values.tolist()

        # validate using the ANUGA tools for complex geometry
        if is_complex(polygon, closed, verbose=verbose):
            # error message taken from ANUGA read_polygon file
            msg = "ERROR: Self-intersecting polygon detected in file "
            msg += filename + ". A complex polygon will not "
            msg += "necessarily break the algorithms within ANUGA, but it"
            msg += "usually signifies pathological data. Please fix this file."
            raise Exception(msg)

        return polygon

    file_path = file

    # open the file
    f = pd.read_csv(file_path, **kwargs)

    # check if there are multiple polygons defined
    ids = f[id_col]
    unique_ids = np.unique(ids)
    # loop through all unique ids
    polygons = []
    for uid in unique_ids:
        new_f = _parse_and_drop_and_validate(f, uid)
        polygons.append(new_f)
    if len(unique_ids) > 1:
        # multiple polygons defined
        pass
    elif len(unique_ids) == 1:
        # single, strip the list
        polygons = polygons[0]
    else:
        raise ValueError("No polygons found.")

    if write:
        raise NotImplementedError("Needs implementation and tests.")
        index = file_path.find(".csv")
        new_file_path = file_path[:index] + "_clean" + file_path[index:]
        new_f.to_csv(os.path.join(data_dir, new_file_path), index=False, header=False)
    return polygons


def raster_to_xyz(file):
    """Process a raster to xyz.

    There are many available functions to do this, but this works to get
    things into the ANUGA format, so why not rewrite it?!
    """
    topography_file = os.path.join(data_dir, file)
    src = rasterio.open(topography_file)
    topo = src.read().squeeze()  # 2d array of topo
    topo_long = topo.flatten()
    l, b, r, t = src.bounds  # bounding box of image
    res = src.res  # resolution of image
    meshX, meshY = np.meshgrid(
        np.arange(l, r, res[0]), np.arange(t, b, -res[0])
    )  # meshgrid of X and Y
    X = np.array(meshX.flatten())  # flatten X and Y
    Y = np.array(meshY.flatten())
    topo_xyz = np.column_stack((X, Y, topo_long))
    return topo_xyz
