
import setuptools
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='hydrafloods',
    version='0.3.5',
    description='HYDrologic Remote sensing Analysis for Floods',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='http://github.com/servir-mekong/hydra-floods',
    packages=setuptools.find_packages(),
    author='Kel Markert',
    author_email='kel.markert@gmail.com',
    license='GNU GPL v3.0',
    zip_safe=False,
    include_package_data=True,
    entry_points={
    'console_scripts': [
        'hydrafloods = hydrafloods.hfcli:main',
    ]},
    install_requires=[
        'simplecmr',
        'earthengine-api',
        'gcsfs',
        'fire',
        'pipetools'
    ],
)
