!/usr/bin/env python

"""
Adjust the CABLE namelist files on the fly, i.e. the user passes a replacement
dictionary of key=value and this script will change and *replace* the existing
namelist file

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (30.07.2018)"
__email__ = "mdekauwe@gmail.com"

import sys
import os
import shutil
import tempfile

def adjust_nml_file(fname, replacements):
    """
    Adjust the params/flags in the CABLE namelise file. Note this writes
    over whatever file it is given!

    Parameters:
    ----------
    replacements : dictionary
        dictionary of replacement values.
    """
    f = open(fname, 'r')
    param_str = f.read()
    f.close()
    new_str = replace_keys(param_str, replacements)
    fd, path = tempfile.mkstemp()
    os.write(fd, str.encode(new_str))
    os.close(fd)
    shutil.copy(path, fname)
    os.remove(path)

def replace_keys(text, replacements_dict):
    """ Function expects to find CABLE namelist file formatted key = value.
    Parameters:
    ----------
    text : string
        input file data.
    replacements_dict : dictionary
        dictionary of replacement values.
    Returns:
    --------
    new_text : string
        input file with replacement values
    """
    lines = text.splitlines()
    for i, row in enumerate(lines):
        # skip blank lines
        if not row.strip():
            continue
        if "=" not in row:
            lines[i] = row
            continue
        elif not row.startswith("&"):
            key = row.split("=")[0]
            val = row.split("=")[1]
            lines[i] = " ".join((key.rstrip(), "=",
                                 replacements_dict.get(key.strip(),
                                 val.lstrip())))

    return '\n'.join(lines) + '\n'
