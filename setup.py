from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="devkidcc",
    version="0.5.0",
    author="Sean Wilson",
    author_email="sean.wilson@mcri.edu.au",
    description="Python wrapper for DevKidCC: Developing Kidney Cell Classifier",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KidneyRegeneration/DevKidCC-python",
    packages=find_packages(),
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
    install_requires=[
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "scanpy>=1.9.0",
        "anndata>=0.8.0",
        "rpy2>=3.5.0",
        "anndata2ri>=1.1.0",
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