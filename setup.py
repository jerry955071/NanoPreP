from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="nanoprep-ccc",
    version="0.0.3",
    author="Chia-Chen Chu",
    author_email="jerry955071@gmail.com",
    description="A fully-equipped, fast, and memory-efficient pre-processor for ONT transcriptomic data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jerry955071/NanoPreP",
    packages=find_packages(),
    install_requires=[
        "edlib>=1.3.8",
        "numpy"
        ],
    entry_points={
        'console_scripts': [
            'nanoprep = NanoPreP.__main__:main',
        ]
    }
)
