from setuptools import setup, find_packages

setup(name='aacrgenie',
      version='1.7.0-veoibddevelop',
      description='Processing and validation for GENIE',
      url='https://github.com/Sage-Bionetworks/Genie',
      author='Thomas Yu',
      author_email='thomasyu888@gmail.com',
      license='MIT',
      packages=find_packages(),
      zip_safe=False,
      python_requires='>=3.5',
      entry_points={
        'console_scripts': ['genie = genie.__main__:main']},
      scripts=['bin/input_to_database.py'],
      install_requires=[
        'pandas>=0.20.0',
        'synapseclient>=1.9',
        'httplib2>=0.11.3',
        'pycrypto>=2.6.1',
        'PyYAML>=5.1'])
