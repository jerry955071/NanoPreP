from setuptools import setup, find_packages
from NanoPrePro._version import __version__

# read README.md
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="nanoprepro",
    version=__version__,
    author="Chia-Chen Chu",
    author_email="b05b01002@ntu.edu.tw",
    description="A fully-equipped, fast, and memory-efficient pre-processor for ONT transcriptomic data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Woodformation1136/NanoPreP",
    packages=find_packages(),
    install_requires=["edlib>=1.3.8", "numpy", "plotly", "pandas"],
    entry_points={
        "console_scripts": [
            "nanoprepro = NanoPrePro.__main__:main"
        ]
    }
)
