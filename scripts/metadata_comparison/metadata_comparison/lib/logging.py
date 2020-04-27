#!/usr/bin/env python3

import logging

def set_log_verbosity(verbose: bool) -> None:
    if verbose:
        logging.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=logging.INFO)
    else:
        logging.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=logging.WARNING)


def quieten_chatty_imports() -> None:
    logging.getLogger('googleapiclient.discovery_cache').setLevel(logging.ERROR)
    logging.getLogger('googleapiclient.discovery').setLevel(logging.WARNING)
