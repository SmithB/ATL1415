import os
from setuptools import setup, find_packages

# get long_description from README.md
with open("README.md", "r") as fh:
    long_description = fh.read()


# list of all scripts to be included with package
scripts = [os.path.join('scripts',f) for f in os.listdir('scripts') if not (f[0]=='.' or f[-1]=='~' or os.path.isdir(os.path.join('scripts', f)))] +\
    [os.path.join('ATL1415', f) for f in ['ATL11_to_ATL15.py']]

setup(
    name='ATL1415',
    version='1.0.0.0',
    description='Driver for elevation-change mapping in Python.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/bsmith/ATL1415',
    author='Ben Smith',
    author_email='besmith@uw.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
    ],
    keywords='altimetry, least-squares, yak_shaving',
    packages=find_packages(),
    scripts=scripts,
    include_package_data=True,
    package_data={'ATL1415':['resources/*', 'resources/templates/*']}
)
