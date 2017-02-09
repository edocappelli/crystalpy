
__authors__ = ["E Cappelli, M Glass, M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "23/11/2016"

from setuptools import setup

setup(name='crystalpy',
      version='0.0.1',
      description='Python crystal polarization calcution',
      author='Edoardo Cappelli, Mark Glass, Manuel Sanchez del Rio',
      author_email='srio@esrf.eu',
      url='https://github.com/edocappelli/crystalpy/',
      packages=['crystalpy',
                'crystalpy.util',
                'crystalpy.diffraction',
                'crystalpy.polarization',
                'crystalpy.examples',
                'crystalpy.tests'],
      install_requires=[
                        'numpy',
                        'scipy'
                       ],
      # test_suite='tests'
      )
