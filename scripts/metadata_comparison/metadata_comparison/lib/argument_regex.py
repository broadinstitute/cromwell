#!/usr/bin/env python3

import argparse
import re


def workflow_regex_validator(value: str) -> str:
    """Makes sure that a value is a valid Cromwell workflow ID then returns the workflow ID"""
    workflow_regex=re.compile('^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$')
    if not workflow_regex.match(value):
        msg = f'Invalid workflow ID {value}. Expected {workflow_regex.pattern}'
        raise argparse.ArgumentTypeError(msg)
    else:
        return value


def url_regex_validator(value: str) -> str:
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


def gcs_path_regex_validator(value: str) -> (str, str):
    """
    Validates then extracts the bucket and object-path from a GS string. Returned as a pair.
    eg:
        'gs://bucket/path/to/directory/' -> ('bucket', 'path/to/directory')
        or
        'gs://bucket/path/to/file.ext' -> ('bucket', 'path/to/file.ext')
    """
    bucket_class = 'a-zA-Z0-9-'
    object_class = '_\\.' + bucket_class
    gcs_regex = re.compile(f'^gs://(?P<bucket>[{bucket_class}]+)/(?P<object>([{object_class}]+/)*[{object_class}]+)/?$')
    m = gcs_regex.match(value)
    if m:
        return m.group('bucket'), m.group('object')
    else:
        msg = f'Invalid GCS path {value}. Expected {gcs_regex.pattern}'
        raise argparse.ArgumentTypeError(msg)


def digester_version_regex_validator(value: str) -> str:
    """
    Validates that digester version looks like 0.0.1
    """
    digester_version_regex = re.compile('^\\d+\\.\\d+\\.\\d+$')
    m = digester_version_regex.match(value)
    if m:
        return m.group(0)
    else:
        msg = f'Invalid digester version {value}. Expected {digester_version_regex.pattern}'
        raise argparse.ArgumentTypeError(msg)
