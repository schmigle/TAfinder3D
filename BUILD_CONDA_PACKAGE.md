# Building and Distributing TAfinder3D via Conda

## Prerequisites

### Add channels to your conda config (one-time setup):
```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels defaults
```

### Or set channels for this build only:
```bash
export CONDA_CHANNEL_PRIORITY=flexible
```

## Quick Start - Build and Upload to Your Own Channel

### 1. Install conda-build
```bash
conda install -c conda-forge conda-build anaconda-client
```

### 2. Build the package (use ONE of these methods)

**Method A: With channel flags (recommended)**
```bash
cd /Users/moshe.steyn/Downloads/TAfinder3D
conda build . --channel conda-forge --channel bioconda --output-folder ./conda-bld
```

**Method B: If channels already configured**
```bash
cd /Users/moshe.steyn/Downloads/TAfinder3D
conda build . --output-folder ./conda-bld
```

### 3. Upload to Anaconda Cloud
```bash
# Login to anaconda.org
anaconda login

# Upload the built package
anaconda upload ./conda-bld/noarch/tafinder3d-1.0.0-py_0.tar.bz2

# Make it public
anaconda show schmigle/tafinder3d
```

### 4. Users can now install with:
```bash
conda install -c schmigle -c conda-forge -c bioconda tafinder3d
```

## Alternative: Test Locally First

### Install from local build:
```bash
conda install --use-local tafinder3d
```

### Test it works:
```bash
tafinder -h
```

## Troubleshooting

If you get "package does not exist" errors:
1. Make sure channels are added: `conda config --show channels`
2. Use explicit channel flags: `--channel conda-forge --channel bioconda`
3. Check your environment has the channels: `conda config --show-sources`

## Files Created

- **setup.py**: Python package setup
- **meta.yaml**: Conda recipe with all dependencies
- **conda_build_config.yaml**: Build configuration
- **BUILD_CONDA_PACKAGE.md**: This file

## Notes

- Package includes ALL conda dependencies from tafinder3d.yml
- Package includes tafinder_db directory with all search databases (~12MB)
- GPU-accelerated MMseqs2 and Foldseek
- Entry point creates 'tafinder' command automatically
- **PyPI dependencies**: After installing via conda, users must install:
  ```bash
  pip install torch==2.9.0 transformers==4.57.1 sentencepiece==0.2.1 protobuf==6.33.0 bio==1.8.1
  ```
