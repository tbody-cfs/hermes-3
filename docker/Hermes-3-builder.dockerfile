# Build as "hermes-3-builder"
# with sudo docker build -f docker/Hermes-3-builder.dockerfile -t hermes-3-builder docker

# Use a spack image with a pinned SHA
FROM spack/ubuntu-jammy:develop

# Install OS packages needed to build the software
RUN apt-get -yqq update && apt-get -yqq upgrade \
 && apt-get -yqq install --no-install-recommends git build-essential vim cmake \
 && rm -rf /var/lib/apt/lists/*

# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir -p /opt/spack-environment
COPY docker/image_ingredients/docker_spack.yaml /opt/spack-environment/spack.yaml

# Install the software
WORKDIR /opt/spack-environment
RUN spack env activate . && spack install --fail-fast && spack gc -y

# Make an 'entrypoint.sh' script which activates the spack environment
RUN spack env activate --sh -d . > activate.sh
