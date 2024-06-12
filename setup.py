#! /usr/bin/env python
from setuptools import setup

setup(
    name="blueswede",
    version="0.0.1",
    author="Texas A&M Surface Processes",
    license="MIT",
    description="Utilities for working with ANUGA",
    long_description=open("README.rst").read(),
    packages=["blueswede"],
    include_package_data=True,
    url="",
    install_requires=["netCDF4", "scipy", "numpy"],
)
