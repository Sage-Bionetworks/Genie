from setuptools import setup, find_packages

setup(name='aacrgenie',
      version='1.6.2',
      description='Processing and validation for GENIE',
      url='https://github.com/Sage-Bionetworks/Genie',
      author='Thomas Yu',
      author_email='thomasyu888@gmail.com',
      license='MIT',
      packages=find_packages(),
      zip_safe=False,
      data_files=[('genie',['genie/addFeatureType.sh','genie/createGTF.sh'])],
      entry_points = {
        'console_scripts': ['genie = genie.__main__:main']},
      install_requires=[
        'pandas>=0.20.0',
        'synapseclient',
        'httplib2',
        'pycrypto'])
