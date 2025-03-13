docker run -it \
  --user root \
  -v /home/jan/Desktop/nektar:/home/nektar/nektar-src \
  -w /home/nektar/nektar-src \
  nektarpp/nektar-dev:latest /bin/bash

mkdir build
cd build
cmake ..
make -j$(nproc)

IncNavierStokesSolver geometry.xml session.xml

mpiexec -np 1 build/solvers/IncNavierStokesSolver/IncNavierStokesSolver coil.xml session.xml -v