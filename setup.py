from setuptools import setup, find_packages

setup(name='aacrgenie',
      version='1.1.0',
      description='Processing and validation for GENIE',
      url='https://github.com/Sage-Bionetworks/Genie',
      author='Thomas Yu',
      author_email='thomasyu888@gmail.com',
      license='MIT',
      packages=find_packages(),
      zip_safe=False,
      entry_points = {
        'console_scripts': ['genie = processing.__main__:main']},
      install_requires=[
        'pandas>=0.20.0',
        'synapseclient',
        'httplib2',
        'pycrypto'])
