#!/usr/bin/env python

"""
Run CABLE spatially.

This script sets various things within a qsub script and then submits the run.

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (27.07.2018)"
__email__ = "mdekauwe@gmail.com"

import subprocess
import sys

start_yr = 1950
end_yr = 1951
qs_cmd = 'qsub -v start_yr=%d,end_yr=%d %s' % (start_yr, end_yr, template_fn)

error = subprocess.call(qs_cmd, shell=True)
if error is 1:
    raise("Job failed to submit")
