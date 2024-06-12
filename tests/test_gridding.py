# unit tests for gridding

import numpy as np

import pytest
from pathlib import Path

import unittest.mock as mock

import os
import shutil

from blueswede import gridding


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
