from . import example_filetype_format

BASE_CLASS = example_filetype_format.FileTypeFormat

def make_format_registry_dict(cls_list):
    """Use an object's _fileType attribute to make a class lookup dictionary.
    
    Args:
        cls_list: A list of Python classes.
    Returns:
        A dictionary mapping the class._fileType to the class.
    """

    return {cls._fileType: cls for cls in cls_list}

def get_subclasses(cls):
    for subclass in cls.__subclasses__():
        yield from get_subclasses(subclass)
        yield subclass

PROCESS_FILES_LIST = [x for x in get_subclasses(BASE_CLASS)]

PROCESS_FILES = make_format_registry_dict(cls_list=PROCESS_FILES_LIST)
