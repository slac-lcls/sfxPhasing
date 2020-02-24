import os
import imp
from setuptools import setup
from setuptools import find_packages

VERSION = imp.load_source(
        'sfxPhasing.version', os.path.join('sfxPhasing', 'version.py')).VERSION

with open('README.md') as f:
    long_description = f.read()

setup(name='sfxPhasing',
      version=VERSION,
      description='Serial Femtosecond Crystallography Phasing Pipeline.',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='https://github.com/slac-lcls/sfxPhasing',
      author='Chen Wu, Alex Batyuk, Chuck Yoon',
      author_email='yoon82@stanford.edu',
      license='GNU General Public License v2.0',
      packages=find_packages(),
      install_requires=[
          'python'>=2.7.17,
          'numpy'>=1.16.5,
          'pandas'>=0.24.2,
          'seaborn'>=0.9.0,
          'matplotlib',
          'pytest'
      ],
      extras_require={
          'tests': ['pytest'],
          'visualization': ['matplotlib', 'seaborn'],
      },
      classifiers=[
          "Programming Language :: Python :: 2",
          "Operating System :: OS Independent",
      ],
      zip_safe=False)
