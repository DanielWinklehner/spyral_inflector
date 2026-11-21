from setuptools import setup

setup(name='spyral_inflector',
      version='2.3.0',
      description='This module calculates the geometry of a spiral inflector for given initial parameters.',
      url='https://github.com/DanielWinklehner/spyral_inflector',
      author='Daniela Campo, Philip Weigel, Daniel Winklehner',
      author_email='winklehn@mit.edu',
      license='MIT',
      packages=['spyral_inflector', 'spyral_inflector.tracking', 'spyral_inflector.workflows'],
      python_requires='>=3.11',
      entry_points={'console_scripts': ['spyral-inflector = spyral_inflector.workflows.cli:main']},
      zip_safe=False)
