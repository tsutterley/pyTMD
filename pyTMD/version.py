#!/usr/bin/env python
u"""
version.py (11/2023)
Gets semantic version number and commit hash from setuptools-scm
"""
import importlib.metadata

# package metadata
metadata = importlib.metadata.metadata("pyTMD")
# get version
version = metadata["version"]
# append "v" before the version
full_version = "v{0}".format(version)
# get project name
project_name = metadata["Name"]
