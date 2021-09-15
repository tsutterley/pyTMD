FROM python:3.9-slim-buster

LABEL Tyler Sutterley "tsutterl@uw.edu"

ENV DEBIAN_FRONTEND="noninteractive" TZ="America/Los_Angeles"

RUN useradd --create-home --shell /bin/bash pytmd

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        build-essential \
        ca-certificates \
        cmake \
        git \
        libboost-dev \
        libexpat1-dev \
        libldap2-dev \
        libtool \
        liburiparser-dev \
        libxml2-dev  \
        pkg-config \
        wget \
        unzip && \
    apt-get clean

WORKDIR /tmp

ENV JOBS 2

ENV CFLAGS="-fPIC"
ENV ZLIB_VERSION=1.2.11
RUN wget -q http://zlib.net/zlib-${ZLIB_VERSION}.tar.gz && \
    tar -xzf zlib-${ZLIB_VERSION}.tar.gz && \
    cd zlib-${ZLIB_VERSION} && \
    ./configure --prefix=/usr/local && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean && \
    cd contrib/minizip && \
    autoreconf -fvi && \
    ./configure --quiet --prefix=/usr/local/ && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean

ENV CFLAGS="-fPIC"
ENV SZIP_SHORT_VERSION=2.1.1
ENV SZIP_VERSION=2.1.1
RUN wget -q https://support.hdfgroup.org/ftp/lib-external/szip/${SZIP_SHORT_VERSION}/src/szip-${SZIP_VERSION}.tar.gz && \
    tar -xzf szip-${SZIP_VERSION}.tar.gz && \
    cd szip-${SZIP_VERSION} && \
    ./configure --quiet --prefix=/usr/local && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean

ENV OPENSSL_SHORT_VERSION=1.1.1
ENV OPENSSL_VERSION=1.1.1k
RUN wget -q https://www.openssl.org/source/openssl-${OPENSSL_VERSION}.tar.gz && \
    tar -xzf openssl-${OPENSSL_VERSION}.tar.gz && \
    cd openssl-${OPENSSL_VERSION} && \
    ./config shared --prefix=/usr/local && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean

ENV CURL_VERSION=7.77.0
RUN wget -q https://curl.haxx.se/download/curl-${CURL_VERSION}.tar.gz && \
    tar -xzf curl-${CURL_VERSION}.tar.gz && \
    cd curl-${CURL_VERSION} && \
    ./configure --quiet \
        --enable-versioned-symbols \
        --enable-openssl-auto-load-config \
        --with-openssl \
        --with-zlib=/usr/local \
        --prefix=/usr/local && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean

ENV XZ_VERSION=5.2.5
RUN wget -q https://tukaani.org/xz/xz-${XZ_VERSION}.tar.gz && \
    tar -xzf xz-${XZ_VERSION}.tar.gz && \
    cd xz-${XZ_VERSION} && \
    ./configure --quiet --prefix=/usr/local && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean

ENV OPENMPI_SHORT_VERSION=v4.1
ENV OPENMPI_VERSION=4.1.1
RUN wget -q https://download.open-mpi.org/release/open-mpi/${OPENMPI_SHORT_VERSION}/openmpi-${OPENMPI_VERSION}.tar.gz && \
    tar -xzf openmpi-${OPENMPI_VERSION}.tar.gz && \
    cd openmpi-${OPENMPI_VERSION} && \
    ./configure --quiet --prefix=/usr/local && \
    make all install && \
    make clean

ENV SQLITE_VERSION=3350500
ENV SQLITE_YEAR 2021
RUN wget -q https://sqlite.org/${SQLITE_YEAR}/sqlite-autoconf-${SQLITE_VERSION}.tar.gz && \
    tar -xzf sqlite-autoconf-${SQLITE_VERSION}.tar.gz && \
    cd sqlite-autoconf-${SQLITE_VERSION} && \
    ./configure --quiet --prefix=/usr/local && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean

ENV LIBJPEG_SHORT_VERSION=9d
ENV LIBJPEG_VERSION=v9d
RUN wget -q http://ijg.org/files/jpegsrc.${LIBJPEG_VERSION}.tar.gz && \
    tar -xzf jpegsrc.${LIBJPEG_VERSION}.tar.gz && \
    cd jpeg-${LIBJPEG_SHORT_VERSION} && \
    ./configure --quiet --prefix=/usr/local && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean

ENV ZLIBLIB="/usr/local/lib"
ENV ZLIBINC="/usr/local/include"
ENV CPPFLAGS="-I/usr/local/include"
ENV LDFLAGS="-L/usr/local/lib"
ENV LD_LIBRARY_PATH="${ZLIBLIB}:${LD_LIBRARY_PATH}"
ENV CFLAGS="-Wall -O -funroll-loops -malign-loops=2 -malign-functions=2"
ENV LIBPNG_VERSION=1.6.37
RUN wget -q https://download.sourceforge.net/libpng/libpng-${LIBPNG_VERSION}.tar.gz && \
    tar -xzf libpng-${LIBPNG_VERSION}.tar.gz && \
    cd libpng-${LIBPNG_VERSION} && \
    ./configure --quiet --prefix=/usr/local && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean

ENV LIBTIFF_VERSION=4.2.0
RUN wget -q https://download.osgeo.org/libtiff/tiff-${LIBTIFF_VERSION}.tar.gz && \
    tar -xzf tiff-${LIBTIFF_VERSION}.tar.gz && \
    cd tiff-${LIBTIFF_VERSION} && \
    ./configure --quiet --prefix=/usr/local \
        --with-jpeg-include-dir=/usr/local/include \
        --with-jpeg-lib-dir=/usr/local/lib && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean

ENV LIBKML_VERSION=1.3.0
RUN wget -q https://github.com/libkml/libkml/archive/${LIBKML_VERSION}.tar.gz && \
    tar -xzf ${LIBKML_VERSION}.tar.gz && \
    cd libkml-${LIBKML_VERSION} && \
    mkdir build && \
    cd build && \
    cmake \
        -DZLIB_INCLUDE_DIR=/usr/local/include \
        -DMINIZIP_LIBRARY=/usr/local/lib \
        -DMINIZIP_INCLUDE_DIR=/usr/local/lib \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/usr/local/ .. && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean

ENV CC=mpicc
ENV HDF5_VERSION=1.10.5
RUN wget -q https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-${HDF5_VERSION}.tar.gz && \
    tar -xzf hdf5-${HDF5_VERSION}.tar.gz && \
    cd hdf5-${HDF5_VERSION} && \
    ./configure --quiet \
        --enable-parallel \
        --enable-hl \
        --enable-shared \
        --prefix=/usr/local \
        --with-zlib=/usr/local \
        --with-szlib=/usr/local && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean

ENV CC=mpicc
ENV CPPFLAGS="-I/usr/local/include"
ENV LDFLAGS="-L/usr/local/lib -lhdf5_hl -lhdf5 -lz -lm -lcurl"
ENV NETCDF_VERSION=4.7.4
RUN wget -q ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-c-${NETCDF_VERSION}.tar.gz && \
    tar -xzf netcdf-c-${NETCDF_VERSION}.tar.gz && \
    cd netcdf-c-${NETCDF_VERSION} && \
    ./configure --quiet \
        --enable-shared \
        --enable-parallel-tests \
        --enable-netcdf-4 \
        --enable-nczarr \
        --prefix=/usr/local && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean

ENV GEOS_VERSION=3.9.1
RUN wget -q https://download.osgeo.org/geos/geos-${GEOS_VERSION}.tar.bz2 && \
    tar -xjf geos-${GEOS_VERSION}.tar.bz2 && \
    cd geos-${GEOS_VERSION} && \
    ./configure --quiet --prefix=/usr/local && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean

ENV PROJ_VERSION=7.2.1
ENV PROJ_DATUMGRID_VERSION=1.8
ENV PROJ_NETWORK ON
ENV SQLITE3_CFLAGS="-I/usr/local/include"
ENV SQLITE3_LIBS="-L/usr/local/lib -lsqlite3"
RUN wget -q https://download.osgeo.org/proj/proj-${PROJ_VERSION}.tar.gz && \
    wget -q http://download.osgeo.org/proj/proj-datumgrid-${PROJ_DATUMGRID_VERSION}.zip && \
    tar -xzf proj-${PROJ_VERSION}.tar.gz && \
    unzip proj-datumgrid-${PROJ_DATUMGRID_VERSION}.zip -d proj-${PROJ_VERSION}/data/ && \
    cd proj-${PROJ_VERSION} && \
    mkdir build && \
    cd build && \
    cmake \
        -DSQLITE3_INCLUDE_DIR=/usr/local/include/ \
        -DSQLITE3_LIBRARY=/usr/local/lib/libsqlite3.so \
        -DTIFF_INCLUDE_DIR=/usr/local/include \
        -DTIFF_LIBRARY_RELEASE=/usr/local/lib/libtiff.so \
        -DCURL_INCLUDE_DIR=/usr/local/include/ \
        -DCURL_LIBRARY=/usr/local/lib/libcurl.so \
        -DPYTHON_EXECUTABLE=/usr/local/bin/python3 \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/usr/local/ .. && \
    cmake --build . && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean

ENV CC=mpicc
ENV CPPFLAGS="-I/usr/local/include"
ENV LDFLAGS="-L/usr/local/lib"
ENV LIBKML_CFLAGS="/usr/local/include"
ENV LIBKML_LIB="/usr/local/lib"
ENV HDF5_CFLAGS="--enable-parallel --enable-hl --enable-shared"
ENV HDF5_INCLUDE="/usr/local/include"
ENV HDF5_LIBS="/usr/local/lib"
ENV GDAL_VERSION=3.3.0
RUN wget -q https://download.osgeo.org/gdal/CURRENT/gdal-${GDAL_VERSION}.tar.gz && \
    tar -xzf gdal-${GDAL_VERSION}.tar.gz && \
    cd gdal-${GDAL_VERSION} && \
    ./configure --quiet \
        --disable-debug \
        --disable-static \
        --with-hdf5=/usr/local \
        --with-netcdf=/usr/local \
        --with-curl=/usr/local/bin/curl-config \
        --with-crypto=/usr/local \
        --with-geos=/usr/local/bin/geos-config \
        --with-geotiff \
        --with-hide-internal-symbols=yes \
        --with-libkml=/usr/local \
        --with-libkml-inc=/usr/local/include \
        --with-libkml-lib=/usr/local/lib \
        --with-liblzma=/usr/local \
        --with-libtiff=/usr/local \
        --with-libz=/usr/local \
        --with-jpeg=/usr/local \
        --with-openjpeg \
        --with-png=/usr/local \
        --with-proj=/usr/local \
        --with-sqlite3=/usr/local \
        --with-proj=/usr/local \
        --with-rename-internal-libgeotiff-symbols=yes \
        --with-rename-internal-libtiff-symbols=yes \
        --with-threads=yes \
        --without-hdf4 \
        --without-idb \
        --without-jpeg12 \
        --without-perl \
        --without-python \
        --prefix=/usr/local && \
    make --quiet --jobs=${JOBS} && \
    make --quiet install && \
    make clean

RUN ldconfig

RUN rm -rf /tmp/*

WORKDIR /home/pytmd

ENV MPICC=mpicc
ENV CC=mpicc
ENV HDF5_MPI="ON"

RUN pip3 install --no-cache-dir --no-binary=h5py,cartopy \
        cartopy \
        gdal \
        h5py \
        lxml \
        matplotlib \
        mpi4py \
        netCDF4 \
        numpy \
        pandas \
        pyproj \
        python-dateutil \
        pyyaml \
        scikit-learn \
        scipy \
        setuptools_scm \
        zarr && \
    pip3 install --no-cache-dir --no-deps git+https://github.com/tsutterley/read-ICESat-2.git && \
    pip3 install --no-cache-dir --no-deps git+https://github.com/tsutterley/read-ATM1b-QFIT-binary.git

COPY . .

RUN --mount=source=.git,target=.git,type=bind \
    pip install --no-cache-dir --no-deps .

USER pytmd

CMD ["bash"]
