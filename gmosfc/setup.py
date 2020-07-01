#!/usr/bin/env python
import os
import sys
from setuptools import setup


# REQUIREMENTS
# 1. What are the required dependencies?
with open('requirements.txt') as f:
    install_requires = f.read().splitlines()

setup(name='gmosfc',
      version='0.0.8',
      description="Handy GMOS finding chart generator",
      author='Venu Kalari',
      author_email='vkalari@gemini.edu',
      url='https://github.com/astroquackers/gmosfc/',
      license='MIT',
      package_dir={
            'gmosfc': 'gmosfc'},
      packages=['gmosfc'],
      install_requires=install_requires,
      include_package_data=True,
      classifiers=[
          "Development Status :: 4 - Beta",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
          "Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Astronomy",
          ],
      )

