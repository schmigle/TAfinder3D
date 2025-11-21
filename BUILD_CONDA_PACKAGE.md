# Building and Distributing TAfinder3D via Conda

## Quick Start - Build and Upload to Your Own Channel

### 1. Install conda-build
```bash
conda install conda-build anaconda-client
```

### 2. Build the package
```bash
cd /Users/moshe.steyn/Downloads/TAfinder3D
conda build . --output-folder ./conda-bld
```

### 3. Upload to Anaconda Cloud
```bash
# Login to anaconda.org
anaconda login

# Upload the built package
anaconda upload ./conda-bld/noarch/tafinder3d-1.0.0-*.tar.bz2

# Optional: Make it public
anaconda show schmigle/tafinder3d
```

### 4. Users can now install with:
```bash
conda install -c schmigle tafinder3d
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

## Files Created

- **setup.py**: Python package setup (for pip compatibility)
- **meta.yaml**: Conda recipe (defines all dependencies from tafinder3d.yml)
- **BUILD_CONDA_PACKAGE.md**: This file (instructions)

## Notes

- Package includes ALL dependencies from tafinder3d.yml
- GPU-accelerated MMseqs2 and Foldseek
- PyTorch and transformers for structure prediction
- Entry point creates 'tafinder' command automatically

