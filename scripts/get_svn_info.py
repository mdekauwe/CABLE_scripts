#!/usr/bin/env python

"""
Get the SVN url and revision number

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (01.08.2018)"
__email__ = "mdekauwe@gmail.com"


import os

def get_svn_info(here, there):

    os.chdir(there)
    os.system("svn info > tmp_svn")
    fname = 'tmp_svn'
    fp = open(fname, "r")
    svn = fp.readlines()
    fp.close()
    os.remove(fname)

    url = [i.split(":", 1)[1].strip() for i in svn if i.startswith('URL')]
    rev = [i.split(":", 1)[1].strip() for i in svn if i.startswith('Revision')]
    os.chdir(here)

    return url, rev
