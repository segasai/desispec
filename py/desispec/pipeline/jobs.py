#
# See top-level LICENSE.rst file for Copyright information
#
# -*- coding: utf-8 -*-
"""
desispec.pipeline.jobs
==========================

Functions for working with scheduler jobs
"""

from __future__ import absolute_import, division, print_function

import os
import re
import sys

import numpy as np

from desiutil.log import get_logger



def slurm_batch(scripts, deps):
    """Write a dictionary to a file.

    Args:
        path (str): the output file name.
        input (dict): the data.

    Returns:
        nothing.
    """
    for scr in scripts:
        scom = "sbatch {} {}".format(depstr, scr)
        #print("RUN SCRIPTS: {}".format(scom))
        log.debug(time.asctime())
        log.debug(scom)
        sout = sp.check_output(scom, shell=True, universal_newlines=True)
        log.debug(sout)
        p = sout.split()
        jid = re.sub(r'[^\d]', '', p[3])
        jobids.append(jid)
    with open(path, "w") as f:
        ydump(input, f, default_flow_style=False)
    return
