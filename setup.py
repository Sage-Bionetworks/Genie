from setuptools import setup, find_packages

setup(name='veoibd-data-pipeline',
      version='1.6.2-rc1',
      description='Processing and validation for VEOIBD',
      url='https://github.com/veo-ibd/veoibd-data-pipeline',
      author='Kenneth Daily',
      author_email='kenneth.daily@sagebionetworks.org',
      license='MIT',
      packages=find_packages(),
      zip_safe=False,
      python_requires='>=3.5',
      entry_points={
        'console_scripts': ['veoibd = genie.__main__:main']},
      scripts=['bin/input_to_database.py'],
      install_requires=[
        'pandas>=0.20.0',
        'synapseclient>=1.9',
        'httplib2>=0.11.3',
        'pycrypto>=2.6.1',
        'PyYAML>=5.1'])
