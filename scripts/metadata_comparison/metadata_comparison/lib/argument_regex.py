#!/usr/bin/env python3

import argparse
import re

def workflow_regex_validator(value):
    """Makes sure that a value is a valid Cromwell workflow ID then returns the workflow ID"""
    workflow_regex=re.compile('^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$')
    if not workflow_regex.match(value):
        msg = f'Invalid workflow ID {value}. Expected {workflow_regex.pattern}'
        raise argparse.ArgumentTypeError(msg)
    else:
        return value


def url_regex_validator(value):
    """
    Validates then extract the root of the Cromwell URL from the various URL strings which might be provided.
    Deliberately flexible because it's tedious to remember which script requires which type of format.
    eg:
        'http://localhost' => 'http://localhost'
        'http://localhost:8000' => 'http://localhost:8000'
        'http://localhost:8000/' => 'http://localhost:8000'
        'http://localhost:8000/api/workflows/' => 'http://localhost:8000'
        'http://localhost:8000/custom/prefix/api/workflows/' => 'http://localhost:8000/custom/prefix'
    """
    url_regex = re.compile('(http(s?)://((?!/api).)*[^/])(/(api.*)?)?')
    m = url_regex.match(value)
    if m:
        return m.group(1)
    else:
        msg = f'Invalid Cromwell URL {value}. Expected {url_regex.pattern}'
        raise argparse.ArgumentTypeError(msg)


def gcs_path_regex_validator(value):
    """
    Validates then extracts the bucket and object-path from a GS string. Returned as a pair.
    eg:
        'gs://bucket/path/to/directory/' -> ('bucket', 'path/to/directory')
    """
    gcs_regex = re.compile('^gs://([a-zA-Z0-9-]+)/(([a-zA-Z0-9-]+/)*[a-zA-Z0-9-]+)/?$')
    m = gcs_regex.match(value)
    if m:
        return m.group(1), m.group(2)
    else:
        msg = f'Invalid GCS path {value}. Expected {gcs_regex.pattern}'
        raise argparse.ArgumentTypeError(msg)
