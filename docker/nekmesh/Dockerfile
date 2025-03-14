ARG ENV_IMAGE=nektarpp/nektar-env:latest
FROM $ENV_IMAGE AS build

LABEL maintainer="Nektar++ Development Team <nektar-users@imperial.ac.uk>"

COPY --chown=nektar:nektar . /home/nektar

ARG INSTALL_PREFIX=/usr/local

RUN mkdir build && cd build &&                                  \
    cmake                                                       \
    -C ../cmake/NekMesh.cmake                                   \
    -DNEKTAR_LIB_DIR=lib                                        \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}                    \
    .. &&                                                       \
    make -j$(nproc) &&                                          \
    make -j$(nproc) user-guide-pdf developer-guide-pdf

USER root
RUN cd build && make install && cd python && python3 -m pip install --break-system-packages . && cd ../..
USER nektar

FROM debian:bookworm-slim
RUN DEBIAN_RELEASE=bookworm && \
    echo "deb http://deb.debian.org/debian ${DEBIAN_RELEASE} non-free" > \
        /etc/apt/sources.list.d/debian-non-free.list && \
    apt-get update && \
    apt-get install -y \
    openssh-client nano vim binutils \
    libboost-iostreams1.74 libboost-system1.74 libboost-program-options1.74 \
    libfftw3-double3 libhdf5-103 libhdf5-openmpi-103 liblapack3 liboce-foundation11 \
    liboce-modeling11 liboce-ocaf11 libtet1.5 libtriangle-1.6 python3-numpy libvtk9.1 \
    openmpi-bin libpetsc-real3.18 libptscotch-7.0 libarpack2 libtinyxml2.6.2v5 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
ARG INSTALL_PREFIX=/usr/local
COPY --from=build ${INSTALL_PREFIX} ${INSTALL_PREFIX}

# Set up entrypoint for copying test files.
COPY docker/nektar/docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh && \
    ln -s /usr/local/bin/docker-entrypoint.sh /
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

RUN groupadd nektar && useradd -m -g nektar nektar && rm /etc/ld.so.cache && ldconfig
USER nektar
WORKDIR /home/nektar
CMD ["/bin/bash"]
