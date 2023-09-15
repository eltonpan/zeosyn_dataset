"""Load pre-defined datasets."""
import gzip
import os

__author__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'
__maintainer__ = 'Haoyan Huo'

__all__ = [
    'open_data'
]

FILEDIR = os.path.dirname(os.path.realpath(__file__))


def open_data(file_name, mode='r'):
    """Open a data file for reading. The file could be gzip'ed."""
    file_name = os.path.join(FILEDIR, 'data', file_name)
    if os.path.exists(file_name):
        return open(file_name, mode)  # pylint: disable=consider-using-with

    if os.path.exists(file_name + '.gz'):
        return gzip.open(file_name + '.gz', mode)

    raise FileNotFoundError(f'Cannot find data file {file_name}')
