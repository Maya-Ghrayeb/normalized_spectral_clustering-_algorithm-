from setuptools import setup, Extension

setup(name='myspkmeans',
      version='1.0',
      description='myspkmeans setup.py',
      ext_modules=[Extension('myspkmeans', sources=['spkmeansmodule.c', 'spkmeans.c'])])