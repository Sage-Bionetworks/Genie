"""genie package setup"""
import os
from setuptools import setup, find_packages

# figure out the version
about = {}
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, "genie", "__version__.py")) as f:
    exec(f.read(), about)

# Add readme
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="aacrgenie",
    version=about["__version__"],
    description="Processing and validation for GENIE",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Sage-Bionetworks/Genie",
    author="Thomas Yu",
    author_email="thomas.yu@sagebionetworks.org",
    license="MIT",
    packages=find_packages(),
    zip_safe=False,
    python_requires=">=3.7,<3.10",
    entry_points={"console_scripts": ["genie = genie.__main__:main"]},
    scripts=["bin/input_to_database.py", "bin/database_to_staging.py"],
    install_requires=[
        "pandas>=1.0",
        "synapseclient>=2.5.1",
        "httplib2>=0.11.3",
        "pycryptodome>=3.12.0",
        "PyYAML>=5.1",
        "chardet>=3.0.4",
    ],
)
