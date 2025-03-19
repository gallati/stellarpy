# StellarPy setup

from setuptools import setup, find_packages

setup(
    name="stellarpy",
    version="0.1",
    packages=find_packages(),
    include_package_data=True,
    package_data={"stellarpy": ["hertzsprung-russell-data.csv"]},
    description="Python package for solving the stellar-interior equations.",
    author="Nuno Cerviño Luridiana",
    author_email="ncervino@ucm.es",
    url="https://github.com/gallati/stellar-interior-numerical-model"
)