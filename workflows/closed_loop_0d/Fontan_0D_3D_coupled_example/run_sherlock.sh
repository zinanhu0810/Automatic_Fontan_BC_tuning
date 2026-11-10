#!/bin/bash
#SBATCH --job-name=fontan_0d3d
#SBATCH --partition=amarsden
#SBATCH --output=fontan_0d3d.o%j
#SBATCH --error=fontan_0d3d.e%j
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=8G

set -euo pipefail

: "${SV_MULTIPHYSICS_BIN:?Set SV_MULTIPHYSICS_BIN to the svmultiphysics executable.}"
: "${SVZERO_INTERFACE_LIB:?Set SVZERO_INTERFACE_LIB to libsvzero_interface.so.}"

module --force purge
module load devel math imkl/2019 cuda/10.2.89 openmpi/4.0.3 hdf5/1.10.6 gcc/12.4.0 \
  freetype/2.8.1 libxml2/2.9.4 libxkbcommon/0.9.1 libfabric/1.9.1 ucx/1.8.1 system \
  szip/2.1.1 cmake/3.24.2 libpng/1.2.57 fontconfig/2.12.4 llvm/4.0.0 x11/7.7 zlib/1.2.11

sed "s|__SVZERO_INTERFACE_LIB__|${SVZERO_INTERFACE_LIB}|g" solver.template.xml > solver.xml
srun "${SV_MULTIPHYSICS_BIN}" solver.xml
