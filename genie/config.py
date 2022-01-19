"""Configuration to obtain registry classes"""
import importlib
import logging

from . import example_filetype_format

logger = logging.getLogger(__name__)

# BASE_CLASS = example_filetype_format.FileTypeFormat


def make_format_registry_dict(cls_list):
    """Use an object's _fileType attribute to make a class lookup dictionary.

    Args:
        cls_list: A list of Python classes.

    Returns:
        A dictionary mapping the class._fileType to the class.

    """
    return {cls._fileType: cls for cls in cls_list}


def get_subclasses(cls):
    """Gets subclasses of modules and classes"""
    for subclass in cls.__subclasses__():
        yield from get_subclasses(subclass)
        yield subclass


def find_subclasses(package_names, base_class):
    """Finds subclasses of a specified base class
    from a list of package names.

    Args:
        package_names: A list of Python package names as strings.
        base_class: A base class to use to search for subclasses

    Returns:
        A list of subclasses that extend from the base class

    """
    matching_classes = []
    for package_name in package_names:
        importlib.import_module(package_name)

    for cls in get_subclasses(base_class):
        logger.debug("checking {cls}.".format(cls=cls))
        cls_module_name = cls.__module__
        cls_pkg = cls_module_name.split(".")[0]
        if cls_pkg in package_names:
            matching_classes.append(cls)
    return matching_classes


def collect_format_types(package_names):
    """Finds subclasses of the example_filetype_format.FileTypeFormat from a
    list of package names.

    Args:
        package_names: A list of Python package names as strings.

    Returns:
        A list of classes that are in the named packages and subclasses of
        example_filetype_format.FileTypeFormat.

    """
    file_format_list = find_subclasses(
        package_names, example_filetype_format.FileTypeFormat
    )
    file_format_dict = make_format_registry_dict(file_format_list)
    return file_format_dict


# PROCESS_FILES_LIST = [x for x in get_subclasses(BASE_CLASS)]
# PROCESS_FILES = make_format_registry_dict(cls_list=PROCESS_FILES_LIST)
