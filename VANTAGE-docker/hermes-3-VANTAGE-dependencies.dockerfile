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
# copy the files in the hermes-3 repo needed to install dependencies
COPY ./external /root/hermes-3/external
COPY ./spack.yaml /root/hermes-3/spack.yaml
# set the workdir
WORKDIR /root/hermes-3

# Install hermes-3 dependencies via spack
# Activate the hermes-3 environment, install dependencies
RUN <<EOF 
spack env activate . -v gcc
# Concretize (in the case the external dir above is a local dir
# where an `$ spack install` has already taken place)
spack concretize -f
# Install the dependencies
spack install -j 4 --only dependencies
# Uninstall the top-level packages that we expect to develop regularly
spack uninstall -y --dependents boutpp@develop || true
spack uninstall -y --dependents vantagereactions@working || true
spack uninstall -y --dependents neso-particles@working || true
EOF
# set workdir back to /root
WORKDIR /root
# remove hermes-3 clone, leaving dependencies installed
RUN rm -rf /root/hermes-3

ENTRYPOINT ["/bin/bash"]
