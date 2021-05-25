FROM ubuntu:20.04

LABEL Tyler Sutterley "tsutterl@uw.edu"

ENV DEBIAN_FRONTEND="noninteractive" TZ="America/Los_Angeles"

RUN useradd --create-home --shell /bin/bash pytmd

RUN apt-get update -y && \
    apt-get install -y \
        openmpi-bin \
        mpi-default-bin \
        openmpi-doc \
        libopenmpi-dev \
        python3-dev \
        python3-pip \
        gdal-bin \
        libgdal-dev \
        libproj-dev \
        proj-data \
        proj-bin \
        libgeos-dev \
        libhdf5-openmpi-dev \
        libnetcdf-dev \
        git && \
    apt-get clean

WORKDIR /home/pytmd

ENV MPICC=mpicc
ENV CC=mpicc
ENV HDF5_MPI="ON"

RUN pip3 install --no-cache-dir --no-binary=h5py,cartopy \
        setuptools_scm \
        mpi4py \
        numpy \
        scipy \
        lxml \
        pyproj \
        python-dateutil \
        pyyaml \
        pandas \
        scikit-learn \
        matplotlib \
        gdal \
        netCDF4 \
        zarr \
        h5py \
        cartopy && \
    pip3 install --no-cache-dir --no-deps git+https://github.com/tsutterley/read-ICESat-2.git && \
    pip3 install --no-cache-dir --no-deps git+https://github.com/tsutterley/read-ATM1b-QFIT-binary.git

COPY . .

RUN --mount=source=.git,target=.git,type=bind \
    pip install --no-cache-dir --no-deps .

USER pytmd

CMD ["bash"]
