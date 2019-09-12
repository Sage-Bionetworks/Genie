"""genie package setup"""
from setuptools import setup, find_packages

setup(name='aacrgenie',
      version='3.2',
      description='Processing and validation for GENIE',
      url='https://github.com/Sage-Bionetworks/Genie',
      author='Thomas Yu',
      author_email='thomas.yu@sagebionetworks.org',
      license='MIT',
      packages=find_packages(),
      zip_safe=False,
      entry_points={'console_scripts': ['genie = genie.__main__:main']},
      install_requires=['pandas>=0.20.0',
                        'synapseclient>=1.9',
                        'httplib2>=0.11.3',
                        'pycrypto>=2.6.1',
                        'PyYAML>=3.11'])
