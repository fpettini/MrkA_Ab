#!/usr/bin/env bash
# Compila GROMACS 2024.6 con MPI vero, CUDA 12.8 e patch PLUMED 2.10.2 (serve per -hrex).
# PLUMED, Open MPI, cmake e compilatori stanno in un ambiente conda dedicato. CUDA e' quello di
# sistema se $CUDA_HOME/bin/nvcc esiste, altrimenti il CUDA 12.8 di conda nello stesso ambiente.
# Installa in una cartella tua: niente root, niente torch.
#
# Perche' 2024.6 e non 2024.3: la 2024.3 non compila per le GPU Blackwell RTX (sm_120, errore di ptxas
# "threads per SM ... out of range" nel kernel di pruning), corretto nella 2024.6. La 2024.4 corregge
# anche interazioni mancanti vicino al cutoff con il pruning su GPU (issue 5138).
# PLUMED 2.10.2 ha la patch con -hrex solo per la 2024.3: sulla 2024.6 si applica senza conflitti
# (spostamenti di una riga) e le modifiche di PLUMED risultano identiche.
#
# Uso:  bash rest2/00_compila_gromacs_plumed.sh [cartella_installazione]
#       default /mnt/f.pettini/sw ; NJ=64 per i processi di make,
#       CUDA_HOME=... per un CUDA di sistema diverso da /usr/local/cuda-12.8
#       (non 13.x: GROMACS 2024 non lo supporta)
# Alla fine: <cartella>/gromacs-2024.6-plumed/bin/gmx_mpi, e rest2/attiva_gmx2024.sh per usarlo.
set -eo pipefail   # niente -u: gli script di attivazione di conda e GMXRC usano variabili non definite
PREFIX=${1:-/mnt/f.pettini/sw}
GMXV=2024.6
PATCH_ENGINE=gromacs-2024.3
ENV=gmx2024-plumed
NJ=${NJ:-64}
CUDA_HOME=${CUDA_HOME:-/usr/local/cuda-12.8}
QUI=$(cd "$(dirname "$0")" && pwd)

if [ -x "$CUDA_HOME/bin/nvcc" ]; then
  CUDA_CONDA=0
  "$CUDA_HOME/bin/nvcc" --version | tail -1
else
  echo "In $CUDA_HOME/bin non c'e' nvcc: uso il CUDA 12.8 di conda (ambiente ${ENV}-cuda)"
  CUDA_CONDA=1
  ENV=${ENV}-cuda
fi

# architettura delle GPU (es. 12.0 -> 120): va data esplicitamente, la lista di default non la contiene
SM=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | tr -d ' .' | sort -u | paste -sd';' -)
echo "GPU compute capability: $SM"

# 1. ambiente conda, solo conda-forge. I compilatori sono quelli di conda (gcc 14) perche'
#    PLUMED e Open MPI di conda-forge sono compilati con quelli; nvcc 12.8 li accetta.
source "$(conda info --base)/etc/profile.d/conda.sh"
PKG=("plumed=2.10.2=mpi_openmpi*" "openmpi=5" "cmake=3.31" make perl fftw
     c-compiler cxx-compiler "gcc_linux-64=14" "gxx_linux-64=14")
[ $CUDA_CONDA = 1 ] && PKG+=("cuda-version=12.8" cuda-nvcc cuda-cudart-dev libcufft-dev
                             cuda-profiler-api cuda-nvtx-dev cuda-driver-dev)
if ! conda env list | awk '{print $1}' | grep -qx "$ENV"; then
  conda create -y -n "$ENV" -c conda-forge --override-channels "${PKG[@]}"
fi
conda activate "$ENV"
if [ $CUDA_CONDA = 1 ]; then
  CUDA_HOME=$CONDA_PREFIX
  export CUDA_PATH=$CONDA_PREFIX/targets/x86_64-linux   # serve a FindCUDA per trovare gli header
  "$CUDA_HOME/bin/nvcc" --version | tail -1
elif conda list "^cuda-nvcc$" | grep -q "^cuda-nvcc "; then
  echo "L'ambiente $ENV contiene il CUDA di conda (versione precedente di questo script)."
  echo "Rimuovilo con:  conda env remove -n $ENV -y   e rilancia lo script."
  exit 1
fi
plumed config has mpi || { echo "PLUMED senza MPI: ambiente sbagliato"; exit 1; }

# 2. sorgenti GROMACS e patch PLUMED
mkdir -p "$PREFIX/src"
cd "$PREFIX/src"
if [ ! -d gromacs-$GMXV ]; then
  curl -fL -o gromacs-$GMXV.tar.gz https://ftp.gromacs.org/gromacs/gromacs-$GMXV.tar.gz
  tar xzf gromacs-$GMXV.tar.gz
fi
cd gromacs-$GMXV
if [ ! -f Plumed.cmake ]; then
  plumed patch -p -e $PATCH_ENGINE --shared
  if find . -name "*.rej" | grep -q .; then
    echo "La patch PLUMED non si e' applicata del tutto (file .rej in $PWD): mandami l'output"; exit 1
  fi
fi

# 3. configurazione
rm -rf build && mkdir build && cd build
cmake .. \
  -DCMAKE_INSTALL_PREFIX="$PREFIX/gromacs-$GMXV-plumed" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" \
  -DCMAKE_CUDA_COMPILER="$CUDA_HOME/bin/nvcc" -DCMAKE_CUDA_HOST_COMPILER="$CXX" -DCUDA_HOST_COMPILER="$(command -v "$CXX")" \
  -DCUDA_TOOLKIT_ROOT_DIR="$CUDA_HOME" -DCUDAToolkit_ROOT="$CUDA_HOME" \
  -DGMX_MPI=ON -DGMX_GPU=CUDA -DGMX_CUDA_TARGET_SM="$SM" \
  -DGMX_FFT_LIBRARY=fftw3 -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" -DBUILD_TESTING=OFF 2>&1 | tee cmake.log | grep -E "SIMD|CUDA:|MPI:|Error" || true
[ -f Makefile ] || { echo "cmake fallito: vedi $PWD/cmake.log"; exit 1; }

# 4. compilazione e installazione
make -j"$NJ" > make.log 2>&1 || { tail -20 make.log; echo "make fallito: vedi $PWD/make.log"; exit 1; }
make install > install.log 2>&1

# 5. file per attivare l'ambiente quando serve
cat > "$QUI/attiva_gmx2024.sh" <<EOT
# uso: source rest2/attiva_gmx2024.sh
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate $ENV
source "$PREFIX/gromacs-$GMXV-plumed/bin/GMXRC"
EOT
[ $CUDA_CONDA = 0 ] && echo "export LD_LIBRARY_PATH=$CUDA_HOME/lib64\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}" >> "$QUI/attiva_gmx2024.sh"
source "$QUI/attiva_gmx2024.sh"
gmx_mpi --version 2>/dev/null | grep -E "GROMACS version|MPI library|GPU support|SIMD|CUDA runtime|CUDA compiler"
gmx_mpi mdrun -h 2>/dev/null | grep -q "\[no\]hrex" && echo "OK: gmx_mpi mdrun ha -hrex" || { echo "ERRORE: -hrex assente"; exit 1; }
