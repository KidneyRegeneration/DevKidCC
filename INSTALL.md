# DevKidCC Python Wrapper - Installation Guide

Complete installation guide for the DevKidCC Python wrapper (rpy2 approach).

## Quick Installation (Recommended)

```bash
# 1. Install R (if not already installed)
# See platform-specific instructions below

# 2. Install Python package
pip install devkidcc

# 3. Done! R packages will auto-install on first use
```

## Platform-Specific Instructions

### macOS

#### Install R
```bash
# Using Homebrew (recommended)
brew install r

# Or download from CRAN
# https://cran.r-project.org/bin/macosx/
```

#### Install Python Package
```bash
# Create virtual environment (recommended)
python3 -m venv devkidcc-env
source devkidcc-env/bin/activate

# Install devkidcc
pip install devkidcc

# Or from source
git clone https://github.com/KidneyRegeneration/DevKidCC-python
cd DevKidCC-python
pip install -e .
```

#### Common macOS Issues

**Issue: rpy2 won't compile**
```bash
# Install pkg-config
brew install pkg-config

# Set environment variable
export PKG_CONFIG_PATH="/usr/local/lib/pkgconfig:/opt/homebrew/lib/pkgconfig"

# Try installing rpy2 again
pip install rpy2
```

**Issue: R not found**
```bash
# Add R to PATH
echo 'export PATH="/usr/local/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc

# Verify R is accessible
R --version
```

### Linux (Ubuntu/Debian)

#### Install R
```bash
# Update package list
sudo apt-get update

# Install R and development packages
sudo apt-get install r-base r-base-dev

# Install additional dependencies
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev
```

#### Install Python Package
```bash
# Create virtual environment
python3 -m venv devkidcc-env
source devkidcc-env/bin/activate

# Install devkidcc
pip install devkidcc
```

#### Common Linux Issues

**Issue: rpy2 compilation fails**
```bash
# Install R development headers
sudo apt-get install r-base-dev

# Install Python development headers
sudo apt-get install python3-dev

# Try again
pip install rpy2
```

### Windows

#### Install R
1. Download R from https://cran.r-project.org/bin/windows/base/
2. Run the installer
3. Add R to PATH:
   - Open "Environment Variables"
   - Add `C:\Program Files\R\R-4.x.x\bin` to PATH

#### Install Rtools (required for package compilation)
1. Download from https://cran.r-project.org/bin/windows/Rtools/
2. Install to default location

#### Install Python Package
```bash
# Create virtual environment
python -m venv devkidcc-env
devkidcc-env\Scripts\activate

# Install devkidcc
pip install devkidcc
```

#### Common Windows Issues

**Issue: R not found**
- Verify R is in PATH: `R --version` in Command Prompt
- Add R bin directory to PATH manually

**Issue: Rtools not found**
- Ensure Rtools is installed
- Add Rtools to PATH: `C:\rtools43\usr\bin`

## Verifying Installation

### Check Python Package
```python
python -c "import devkidcc; print('DevKidCC wrapper installed!')"
```

### Check R Installation
```bash
R --version
```

### Check R Packages (optional)
```bash
# The wrapper will install these automatically, but you can check:
R -e "library(Seurat)"
R -e "library(DevKidCC)"
```

### Full Test
```python
import scanpy as sc
import numpy as np
from devkidcc import DevKidCCClassifier

# Create dummy data
adata = sc.AnnData(
    X=np.random.randn(100, 2000),
    var=pd.DataFrame(index=[f"Gene_{i}" for i in range(2000)])
)

# This will check everything:
# 1. R is accessible
# 2. R packages can be installed
# 3. Conversion works
# 4. Classification runs
try:
    classifier = DevKidCCClassifier(verbose=True, install_deps=True)
    print("✓ Installation successful!")
except Exception as e:
    print(f"✗ Installation issue: {e}")
```

## Manual R Package Installation

If automatic installation fails, install R packages manually:

```r
# Open R console
R

# Install from CRAN
install.packages("Seurat")
install.packages("devtools")
install.packages("hdf5r")

# Install from GitHub
devtools::install_github("mojaveazure/seurat-disk")
devtools::install_github("powellgenomicslab/scPred")
devtools::install_github("KidneyRegeneration/DevKidCC")
```

Then use Python wrapper with `install_deps=False`:
```python
classifier = DevKidCCClassifier(install_deps=False)
```

## Installing from Source

```bash
# Clone repository
git clone https://github.com/KidneyRegeneration/DevKidCC-python
cd DevKidCC-python

# Install in development mode
pip install -e .

# Install development dependencies
pip install -e ".[dev]"
```

## Docker Installation (Advanced)

For reproducible environments:

```dockerfile
# Dockerfile
FROM rocker/r-ver:4.3.0

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

# Install R packages
RUN R -e "install.packages(c('Seurat', 'devtools', 'hdf5r'))"
RUN R -e "devtools::install_github('mojaveazure/seurat-disk')"
RUN R -e "devtools::install_github('KidneyRegeneration/DevKidCC')"

# Install Python packages
RUN pip3 install devkidcc

# Set working directory
WORKDIR /workspace

CMD ["python3"]
```

Build and run:
```bash
docker build -t devkidcc .
docker run -it -v $(pwd):/workspace devkidcc
```

## Conda Installation (Alternative)

```bash
# Create conda environment
conda create -n devkidcc python=3.10 r-base=4.3
conda activate devkidcc

# Install R packages via conda (some available)
conda install -c conda-forge r-seurat

# Install remaining R packages in R
R -e "devtools::install_github('KidneyRegeneration/DevKidCC')"

# Install Python package
pip install devkidcc
```

## Upgrading

```bash
# Upgrade Python package
pip install --upgrade devkidcc

# Update R packages (in R)
R -e "devtools::install_github('KidneyRegeneration/DevKidCC', force=TRUE)"
```

## Uninstalling

```bash
# Remove Python package
pip uninstall devkidcc

# Remove R packages (optional, in R)
R -e "remove.packages(c('DevKidCC', 'SeuratDisk', 'Seurat'))"
```

## Troubleshooting

### General Debugging

1. **Check R accessibility**:
   ```bash
   which R
   R --version
   ```

2. **Check Python can find R**:
   ```python
   import rpy2.robjects as ro
   print(ro.r('R.version.string'))
   ```

3. **Check rpy2 installation**:
   ```python
   import rpy2
   print(rpy2.__version__)
   ```

4. **Enable verbose output**:
   ```python
   classifier = DevKidCCClassifier(verbose=True)
   ```

### Getting Help

1. Check error messages carefully
2. Search [GitHub Issues](https://github.com/KidneyRegeneration/DevKidCC-python/issues)
3. Create new issue with:
   - Operating system
   - Python version: `python --version`
   - R version: `R --version`
   - Error message
   - Minimal reproducible example

## System Requirements

**Minimum:**
- Python 3.8+
- R 4.0+
- 4GB RAM
- 2GB disk space

**Recommended:**
- Python 3.10+
- R 4.3+
- 8GB+ RAM
- 5GB disk space (for R package cache)

## Next Steps

After successful installation:

1. Try the quick start example in README.md
2. Run example scripts in `examples/`
3. Read the full documentation
4. Classify your own data!

## Support

- Documentation: [README.md](README.md)
- Examples: [examples/](examples/)
- Issues: https://github.com/KidneyRegeneration/DevKidCC-python/issues
- Email: sean.wilson@mcri.edu.au