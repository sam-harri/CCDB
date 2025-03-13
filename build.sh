docker run -it --user root \
  -v /home/jan/Desktop/nektar:/home/nektar/nektar-src \
  -w /home/nektar/nektar-src \
  nektarpp/nektar-dev:latest /bin/bash -c "
  
  cd build && make -j\$(nproc) && cd ..

  apt update && apt install -y curl

  su -c 'curl -LsSf https://astral.sh/uv/install.sh | sh' nektar

  su -c '/home/nektar/.local/bin/uv run sim.py' nektar
"
