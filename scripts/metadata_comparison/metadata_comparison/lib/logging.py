#!/usr/bin/env python3

import logging as log

def set_log_verbosity(verbose):
    if verbose:
        log.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=log.INFO)
    else:
        log.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=log.WARNING)
