# Build as "hermes-3"
# with sudo docker build -f docker/Hermes-3.dockerfile -t hermes-3 docker

# Load the hermes-3-builder image
FROM ghcr.io/boutproject/hermes-3-builder AS builder

# Make a bare image for building hermes-3
FROM ubuntu:22.04

COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software

# paths.view is a symlink, so copy the parent to avoid dereferencing and duplicating it
COPY --from=builder /opt/views /opt/views

# Install useful tools you'll want in the container
RUN apt-get -yqq update && apt-get -yqq upgrade \
 && apt-get -yqq install --no-install-recommends git build-essential vim cmake \
 && rm -rf /var/lib/apt/lists/*

# Change into the /hermes_project, and define paths
WORKDIR /hermes_project
# Build dirs hold the actual executables
ENV HERMES_BUILD_DIR=/hermes_project/build/hermes-3-build
ENV BOUTPP_BUILD_DIR=/hermes_project/build/boutpp-build
ENV HERMES_BUILD_DIR_OVERRIDE=/hermes_project/work/hermes-3-build
ENV BOUTPP_BUILD_DIR_OVERRIDE=/hermes_project/work/boutpp-build
# Source dirs hold the source code. If you want to edit the source,
# do it in the OVERRIDE directories
ENV HERMES_SRC_DIR=/hermes_project/src/hermes-3
ENV BOUTPP_SRC_DIR=/hermes_project/src/hermes-3/external/BOUT-dev
ENV HERMES_SRC_DIR_OVERRIDE=/hermes_project/work/hermes-3
ENV BOUTPP_SRC_DIR_OVERRIDE=/hermes_project/work/BOUT-dev
# Config files set the options used by the CMake build. If you want to
# change the options, do it in the OVERRIDE files
ENV HERMES_CONFIG=/hermes_project/config/hermes_config.cmake
ENV BOUTPP_CONFIG=/hermes_project/config/boutpp_config.cmake
ENV HERMES_CONFIG_OVERRIDE=/hermes_project/work/hermes_config.cmake
ENV BOUTPP_CONFIG_OVERRIDE=/hermes_project/work/boutpp_config.cmake

# Copy in required files for a minimal build of Hermes-3 and BOUT++
COPY . ${HERMES_SRC_DIR}
# Initialize the git submodules (needed for CI/CD build)
RUN git -C ${HERMES_SRC_DIR} submodule update --init --recursive

COPY docker/image_ingredients/enable_c.patch ${BOUTPP_SRC_DIR}/enable_c.patch
RUN git -C ${BOUTPP_SRC_DIR} apply ./enable_c.patch

# Copy in the default CMake config files
COPY docker/image_ingredients/boutpp_config.cmake ${BOUTPP_CONFIG}
COPY docker/image_ingredients/hermes_config.cmake ${HERMES_CONFIG}

# Configure and build BOUT++
RUN . /opt/spack-environment/activate.sh \
&&  cmake -B ${BOUTPP_BUILD_DIR} \
          -S ${BOUTPP_SRC_DIR} \
          -C ${BOUTPP_CONFIG} \
          -Wno-dev \
&& cmake --build ${BOUTPP_BUILD_DIR} --parallel

# Configure and build Hermes
RUN . /opt/spack-environment/activate.sh \
&&  cmake -B ${HERMES_BUILD_DIR} \
          -S ${HERMES_SRC_DIR} \
          -C ${HERMES_CONFIG} \
          -DCMAKE_PREFIX_PATH=${BOUTPP_BUILD_DIR} \
          -Wno-dev \
&& cmake --build ${HERMES_BUILD_DIR} --parallel

# Copy in some helpful commands which can be used in
# the image. Make sure these can be executed when setting
# --user=${UID}:${GID}
COPY docker/image_ingredients/docker_image_commands.sh /bin/image
RUN chmod u+x,o+x /bin/image

# Prevent an issue where only one process can open a HDF file at a time
ENV HDF5_USE_FILE_LOCKING=False

# Copy in a script which runs before any instance is
# launched
COPY docker/image_ingredients/docker_entrypoint.sh /entrypoint.sh
RUN chmod a+x /entrypoint.sh

# Set the default entrypoint and command for the image
ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash"]
