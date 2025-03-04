# unit tests for gridding

import numpy as np

import pytest
from pathlib import Path

import unittest.mock as mock

import os
import shutil

from blueswede import gridding

import netCDF4 as nc


# sample_data_path = os.path.join("test_data", "channel3.sww")


@pytest.fixture(scope="module")
def test_path(request):
    """Return the directory of the currently running test script"""

    return request.path.parent


class Test_files_folders:
    def test_output_none(self, tmp_path: Path, test_path) -> None:
        """Test from the sample data, processing into a dataset."""
        # set up the paths and copy the test data to the temp folder
        sample_data_path = test_path.joinpath("test_data", "channel3.sww")
        _output_path = os.path.join(tmp_path)  # where we will look for the file
        shutil.copy(sample_data_path, _output_path)
        new_sample_data_path = os.path.join(_output_path, "channel3.sww")
        assert os.path.exists(os.path.join(_output_path, "channel3.sww"))

        # run the gridding
        gridding.grid_sww_to_netcdf(sww_file=new_sample_data_path)

        # check the output
        assert os.path.exists(os.path.join(_output_path, "channel3.nc"))

    def test_output_folder(self, tmp_path: Path, test_path) -> None:
        """Test from the sample data, processing into a dataset."""
        # set up the paths
        sample_data_path = test_path.joinpath("test_data", "channel3.sww")
        _output_path = os.path.join(tmp_path)  # just pass the temp folder, no filename

        # run the gridding
        gridding.grid_sww_to_netcdf(sww_file=sample_data_path, nc_file=_output_path)

        # check the output
        assert os.path.exists(os.path.join(_output_path, "channel3.nc"))

    def test_output_folder_filename(self, tmp_path: Path, test_path) -> None:
        """Test from the sample data, processing into a dataset."""
        # set up the paths
        sample_data_path = test_path.joinpath("test_data", "channel3.sww")
        _output_path = os.path.join(tmp_path, "adifferentfilename.nc")

        # run the gridding
        gridding.grid_sww_to_netcdf(sww_file=sample_data_path, nc_file=_output_path)

        # check the output
        assert os.path.exists(os.path.join(tmp_path, "adifferentfilename.nc"))
        assert not os.path.exists(os.path.join(_output_path, "adifferentfilename.nc"))


class Test_nc_formatting:
    def test_utm(self):
        """
        Check the output of a netcdf file that should be UTM referenced.
        """
        pass

    def test_arbitrary(self):
        """
        Check the output of an arbitraty coordinate system.
        """
        # This would be a test for the channel3 data.
        pass

    def test_no_description(self, tmp_path: Path, test_path) -> None:
        """Test from the sample data, processing into a dataset."""
        # set up the paths and copy the test data to the temp folder
        sample_data_path = test_path.joinpath("test_data", "channel3.sww")
        _output_path = os.path.join(tmp_path)  # where we will look for the file
        shutil.copy(sample_data_path, _output_path)
        new_sample_data_path = os.path.join(_output_path, "channel3.sww")

        # run the gridding
        gridding.grid_sww_to_netcdf(sww_file=new_sample_data_path)

        # check the output
        ds = nc.Dataset(os.path.join(_output_path, "channel3.nc"))
        assert not hasattr(ds, "description")

    def test_description(self, tmp_path: Path, test_path) -> None:
        """Test from the sample data, processing into a dataset."""
        # set up the paths and copy the test data to the temp folder
        sample_data_path = test_path.joinpath("test_data", "channel3.sww")
        _output_path = os.path.join(tmp_path)  # where we will look for the file
        shutil.copy(sample_data_path, _output_path)
        new_sample_data_path = os.path.join(_output_path, "channel3.sww")

        # run the gridding
        gridding.grid_sww_to_netcdf(
            sww_file=new_sample_data_path, nc_description="test description"
        )

        # check the output
        ds = nc.Dataset(os.path.join(_output_path, "channel3.nc"))
        assert ds.description == "test description"


class Test_nc_xy_settings:
    def test_default(self, tmp_path: Path, test_path) -> None:
        """Test for no xy spec and grid spacing 10 m"""
        # set up the paths and copy the test data to the temp folder
        sample_data_path = test_path.joinpath("test_data", "channel3.sww")
        _output_path = os.path.join(tmp_path)  # where we will look for the file
        shutil.copy(sample_data_path, _output_path)
        new_sample_data_path = os.path.join(_output_path, "channel3.sww")

        # run the gridding
        gridding.grid_sww_to_netcdf(sww_file=new_sample_data_path)

        # check the output
        ds = nc.Dataset(os.path.join(_output_path, "channel3.nc"))
        assert (
            ds["easting"][1] - ds["easting"][0]
        ) == 10  # check here for the grid spacing

    def test_20m(self, tmp_path: Path, test_path) -> None:
        """Test for no xy spec and grid spacing 20 m"""
        # set up the paths and copy the test data to the temp folder
        sample_data_path = test_path.joinpath("test_data", "channel3.sww")
        _output_path = os.path.join(tmp_path)  # where we will look for the file
        shutil.copy(sample_data_path, _output_path)
        new_sample_data_path = os.path.join(_output_path, "channel3.sww")

        # run the gridding
        gridding.grid_sww_to_netcdf(sww_file=new_sample_data_path, dx=1)

        # check the output
        ds = nc.Dataset(os.path.join(_output_path, "channel3.nc"))
        assert (
            ds["easting"][1] - ds["easting"][0]
        ) == 1  # check here for the grid spacing

    def test_xy_tuple(self, tmp_path: Path, test_path) -> None:
        """Test for given xy spec and grid spacing 1000m

        note: when xy is given, grid spacing should be ignored.
        """
        """Test for no xy spec and grid spacing 20 m"""
        # set up the paths and copy the test data to the temp folder
        sample_data_path = test_path.joinpath("test_data", "channel3.sww")
        _output_path = os.path.join(tmp_path)  # where we will look for the file
        shutil.copy(sample_data_path, _output_path)
        new_sample_data_path = os.path.join(_output_path, "channel3.sww")

        # make the custom grid (choose something weird to see effect)
        xvect, yvect = np.linspace(5, 10, num=23), np.linspace(3, 5, num=13)

        # run the gridding
        gridding.grid_sww_to_netcdf(
            sww_file=new_sample_data_path, dx=20, xy=(xvect, yvect)
        )  # note, the 20 should be ignored

        # check the output
        ds = nc.Dataset(os.path.join(_output_path, "channel3.nc"))
        assert (
            ds["easting"][1] - ds["easting"][0]
        ) != 20  # check here that the input grid spacing was ignored
        assert ds["easting"][0] == 5
        assert ds["northing"][-1] == 5
        assert len(ds["northing"]) == 13


class Test_mask:
    # no sample data for testing
    def test_mask_enabled(self):
        pass
