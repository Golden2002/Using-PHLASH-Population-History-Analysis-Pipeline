# PHLASH Pipeline Installation Guide

This repository provides a fully pre-built Python environment designed for running the PHLASH pipeline and related population genetics analyses.  
The environment is packaged as a compressed virtual environment directory (`phlash_cpu_env_cpu_only.tar.gz`) to ensure complete reproducibility across different Linux systems and HPC clusters.

This distribution is intended for users who wish to run **PHLASH**, **msprime**, **tskit**, **pyslim**, **stdpopsim**, and a number of **JAX-based or GPU-accelerated tools** without the complexity of rebuilding the environment manually.

Two installation options are provided:

1. **Pre-packaged CPU-only virtual environment** (recommended for quick setup)
2. **Manual environment deployment from configuration files** (recommended for flexible or custom installations)

> Pre-built environment packages are available in the [Releases](https://github.com/Golden2002/Using-PHLASH-Population-History-Analysis-Pipeline/releases/tag/v0.0.1-alpha) section of this repository. You can download the latest version from there.
> A provided `dinopy` library is in the reporsitory path `Using-PHLASH-Population-History-Analysis-Pipeline/Installation Guide/dinopy-3.0.0.tar.gz`

---

## Option 1: Pre-packaged CPU-only Virtual Environment

A pre-built CPU-only environment is provided for convenience. This environment contains all necessary packages for PHLASH pipeline execution, excluding GPU dependencies.

### Important Notes:

#### 1. GPU Support Limitations
**The CPU-only environment does NOT include GPU-accelerated dependencies required for running PHLASH on NVIDIA GPUs.** If you need GPU support, you must manually install additional GPU libraries:

```bash
# After activating the environment, install GPU dependencies:
source ~/phlash_cpu_env/bin/activate
pip install cuda-python nvidia-cudnn-cu12 nvidia-cublas-cu12 nvidia-cusolver-cu12 nvidia-cusparse-cu12 nvidia-curand-cu12 nvidia-cufft-cu12 nvidia-cuda-runtime-cu12 nvidia-cuda-nvrtc-cu12 nvidia-cuda-cupti-cu12 nvidia-cudart-cu12 nvidia-cuda-nvcc-cu12 jax[cuda12]
```

#### 2. Dinopy Library Dependency
**The `dinopy` library is not available on PyPI and must be installed separately.** In the pre-built environment, this is already included. However, if you need to reinstall or update it, use the provided distribution:

```bash
# Install dinopy from the provided source distribution
pip install dinopy-3.0.0.tar.gz
```

> **Note:** The `dinopy` package requires Cython for building. If you encounter build errors, install Cython first: `pip install cython`

#### 3. "Which Python" is In the pre-packaged virtual environment

- Python 3.12.5

### Steps

#### 1. Download the environment archive

```bash
wget https://github.com/Golden2002/Using-PHLASH-Population-History-Analysis-Pipeline/releases/tag/v0.0.1-alpha -O ~/phlash_cpu_env_cpu_only.tar.gz
```

#### 2. Extract the environment

```bash
tar -xzf ~/phlash_cpu_env_cpu_only.tar.gz -C ~/
```

#### 3. Activate the environment

```bash
source ~/phlash_cpu_env/bin/activate
```

#### 4. Verify installation

```bash
which python
python -m pip list
```

The environment is now ready to run the PHLASH pipeline and other related analyses.

> **Note:** This CPU-only environment significantly reduces archive size while maintaining all core functionalities. However, for GPU-accelerated computations, additional CUDA libraries must be installed as described above.

---

## Option 2: Manual Installation from Configuration Files

This option allows users to build the environment locally using Python and pip. It is recommended if you need to customize package versions or install on a different system.

### Steps

#### 1. Create and activate a Python virtual environment

```bash
python3 -m venv ~/phlash_cpu_env_manual
source ~/phlash_cpu_env_manual/bin/activate
python -m pip install --upgrade pip setuptools wheel
```

#### 2. Install the `dinopy` package

```bash
# Dinopy is not available on PyPI, install from the provided source distribution
python -m pip install Installation\ Guide/Phlash/dinopy-3.0.0.tar.gz
```

> **Note:** If you encounter build errors during dinopy installation, you may need to install Cython first:
> ```bash
> pip install cython
> ```

#### 3. Install required Python packages

Choose the desired requirements file:

- **Full environment** (includes GPU packages):

  ```bash
  python -m pip install -r /home/litianxing/100My_Jino/114.PHLASH/software/pip_requirements.txt
  ```

- **CPU-only environment** (recommended for initial setup):

  ```bash
  python -m pip install -r /home/litianxing/100My_Jino/114.PHLASH/software/CPU_only_requirements.txt
  ```

#### 4. For GPU support (if using CPU-only environment)

If you selected the CPU-only environment but need GPU acceleration, install the following additional packages:

```bash
# Install NVIDIA CUDA bindings and JAX with CUDA 12 support
pip install cuda-python nvidia-cudnn-cu12 nvidia-cublas-cu12 nvidia-cusolver-cu12 nvidia-cusparse-cu12 nvidia-curand-cu12 nvidia-cufft-cu12 nvidia-cuda-runtime-cu12 nvidia-cuda-nvrtc-cu12 nvidia-cuda-cupti-cu12 nvidia-cudart-cu12 nvidia-cuda-nvcc-cu12 jax[cuda12]

# Verify GPU availability
python -c "import jax; print(jax.devices())"
```

#### 5. Verify installation

```bash
python -m pip list
```

#### 6. Deactivate the environment after use

```bash
deactivate
```

---

## Dinopy Library Installation Details

### About Dinopy
Dinopy is a custom Python library required by the PHLASH pipeline that is not available on PyPI. The source distribution is provided in this repository.

### Installation Methods

#### Method A: From Provided Distribution (Recommended)

```bash
# From the repository root directory
pip install dinopy-3.0.0.tar.gz
```

#### Method B: From Source Directory (if extracted)

```bash
# First extract the source distribution
tar -xzf dinopy-3.0.0.tar.gz
cd dinopy-3.0.0
pip install .
```

### Quick Installation Command

```bash
# One-liner to install dinopy from the provided distribution
curl -L <RELEASE_URL>/dinopy-3.0.0.tar.gz | pip install -
```

### Troubleshooting Dinopy Installation

If you encounter issues installing dinopy:

1. **Cython not installed**:
   ```bash
   pip install cython
   ```

2. **Build dependencies missing**:
   ```bash
   # For Ubuntu/Debian
   sudo apt-get install python3-dev build-essential
   
   # For CentOS/RHEL
   sudo yum install python3-devel gcc gcc-c++
   ```

3. **Permission errors**:
   ```bash
   # Use --user flag or run in virtual environment
   pip install --user Installation\ Guide/Phlash/dinopy-3.0.0.tar.gz
   ```

4. **Python version compatibility**:
   ```bash
   # Check Python version
   python --version
   # Dinopy requires Python 3.7+
   ```

> **Note:** Ensure you are in the repository root directory when using the relative path.

---

## GPU Support Information

### Required GPU Dependencies for PHLASH

To enable GPU acceleration for PHLASH computations, the following dependencies are required:

1. **CUDA Toolkit 12.x** (system-level installation)
2. **NVIDIA CUDA Python bindings**: `cuda-python`
3. **NVIDIA CUDA libraries**: `nvidia-cudnn-cu12`, `nvidia-cublas-cu12`, etc.
4. **JAX with CUDA support**: `jax[cuda12]`

### Installation Notes

- **CUDA Toolkit**: Must be installed system-wide (not via pip). Check with `nvcc --version`.
- **GPU Environment**: The "full environment" option includes these GPU dependencies.
- **CPU-only Environment**: Lacks GPU support; must manually install the above packages.
- **Compatibility**: Ensure CUDA driver version matches the CUDA toolkit version.

### Quick GPU Support Check

After installing GPU dependencies, verify with:

```bash
python -c "
import jax
import jax.numpy as jnp
from jax import random

print('JAX version:', jax.__version__)
print('Available devices:', jax.devices())
print('GPU available:', any(device.platform == 'gpu' for device in jax.devices()))

# Simple GPU test
key = random.PRNGKey(0)
x = random.normal(key, (1000, 1000))
print('GPU test successful:', x.device())
"
```

---

## Notes and Recommendations

- The CPU-only environment is recommended for most users, especially on HPC clusters or systems without GPU support.
- Manual installation is useful for users who need a system-specific Python setup or wish to customize package versions.
- Ensure **Python 3.8 or higher** is used for compatibility with JAX and msprime.
- **Dinopy library** must be installed from the provided source distribution as it is not available on PyPI.
- The pre-built environment includes all necessary dependencies for:
  - **Population genetics simulations and analysis:** msprime, tskit, pyslim, stdpopsim
  - **PHLASH pipeline execution**
  - **Data handling and visualization:** numpy, pandas, matplotlib, plotly
  - **JAX-based computations** (CPU-only or GPU, depending on the selected environment)

### System Requirements

- **CPU-only**: Minimum 8GB RAM, 20GB disk space
- **GPU support**: NVIDIA GPU with CUDA 12.x compatibility, 16GB+ VRAM recommended
- **Operating System**: Linux (Ubuntu 20.04+, CentOS 8+, or compatible distributions)

---

## Troubleshooting

### Common Issues

1. **Dinopy installation fails**: Ensure Cython is installed and build tools are available.
2. **GPU not detected**: Ensure CUDA toolkit is installed system-wide and NVIDIA drivers are up-to-date.
3. **JAX CUDA errors**: Reinstall with `pip install --force-reinstall jax[cuda12]`.
4. **Environment activation issues**: Ensure you're using `bash` or `zsh` shell; for `csh`/`tcsh`, use `source ~/phlash_cpu_env/bin/activate.csh`.

### Getting Help

For installation issues, please:
1. Check the [GitHub Issues](<REPO_URL>/issues) for similar problems
2. Ensure your system meets the minimum requirements
3. Provide error logs and system information when requesting support

---

This guide ensures that users can rapidly deploy the PHLASH pipeline and associated tools in a reproducible environment, minimizing setup errors and maximizing portability across Linux systems and HPC infrastructures.

