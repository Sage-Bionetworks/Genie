# Import logging last to not take in synapseclient logging
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logging.getLogger("keyring").setLevel(logging.WARNING)

# create version in __init__.py
# https://packaging.python.org/en/latest/guides/single-sourcing-package-version/
__version__ = "14.1.2"

__all__ = ["__version__"]
