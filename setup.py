
import setuptools
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='hydrafloods',
      version='0.0.15',
      description='HYDrologic Remote sensing Analysis for Floods',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='http://github.com/servir-mekong/hydra-floods',
      packages=setuptools.find_packages(),
      author='Kel Markert',
      author_email='kel.markert@gmail.com',
      license='MIT',
      zip_safe=False)
