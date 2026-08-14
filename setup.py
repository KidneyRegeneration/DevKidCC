from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="devkidcc",
    version="0.5.1",
    author="Sean Wilson",
    author_email="sean.wilson@mcri.edu.au",
    description="Python wrapper for DevKidCC: Developing Kidney Cell Classifier",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KidneyRegeneration/DevKidCC/tree/python-wrapper",
    packages=find_packages(),
    package_data={
        "devkidcc": ["run_dkcc.R", "data/reference_genes.txt"],
    },
    include_package_data=True,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    # No rpy2 / anndata2ri: the backend talks to R over a subprocess, which is
    # the whole reason this wrapper exists (rpy2 and reticulate segfault when
    # both are loaded). scipy is required — the CSV writer slices a csc_matrix.
    install_requires=[
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        "anndata>=0.8.0",
        "scanpy>=1.9.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.11.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.9",
        ],
    },
)