#!/usr/bin/env python3

import logging
import warnings


def set_log_verbosity(verbose: bool) -> None:
    level = logging.INFO if verbose else logging.WARNING
    logging.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=level)


def quieten_chatty_imports() -> None:
    logging.getLogger('googleapiclient.discovery_cache').setLevel(logging.ERROR)
    logging.getLogger('googleapiclient.discovery').setLevel(logging.WARNING)
    # Controversial and doesn't seem to work for the tests anyway, YMMV.
    # warnings.filterwarnings("ignore", "Your application has authenticated using end user credentials")
