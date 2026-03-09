# Build stage with Spack pre-installed and ready to be used
FROM spack/ubuntu-noble:1.1.0

RUN apt update && \
 apt install -y git --no-install-recommends \
 && rm -rf /var/lib/apt/lists/*
# use develop version of spack repos
RUN sed -i '/^[[:space:]]*branch:/ s|releases/v2025\.11|develop|g' /opt/spack/etc/spack/defaults/base/repos.yaml
# update spack repos
RUN spack repo update
# find the gcc compiler
RUN spack compiler find gcc
# clone the hermes-3 repo
RUN git clone https://github.com/UKAEA-Edge-Code/hermes-3.git hermes-3
# update the submodules
WORKDIR /root/hermes-3
RUN git submodule update --init --recursive

# Install hermes-3 dependencies via spack
# Activate the hermes-3 environment, install dependencies
RUN <<EOF 
spack env activate . -v gcc
# Install the dependencies
spack install -j 4 --only dependencies
EOF
# uninstall any top-level packages
RUN <<EOF 
spack env activate . -v gcc
# Uninstall the packages that we expect to develop regularly
spack uninstall -y --dependents boutpp@develop || true
spack uninstall -y --dependents vantagereactions@working || true
spack uninstall -y --dependents neso-particles@working || true
EOF
# set workdir back to /root
WORKDIR /root
# remove hermes-3 clone, leaving dependencies installed
RUN rm -rf /root/hermes-3

ENTRYPOINT ["/bin/bash"]
