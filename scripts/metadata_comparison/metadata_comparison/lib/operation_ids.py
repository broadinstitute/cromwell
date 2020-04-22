#!/usr/bin/env python3

import re

operation_regex = re.compile('^.*/([0-9]*)')
papi_v1_operation_regex = re.compile('^projects/.*/operations/([^/]*)')
papi_v2alpha1_operation_regex = re.compile('^projects/.*/operations/([^/]*)')
papi_v2alpha1_operation_regex = re.compile('^projects/.*/operations/([^/]*)')



def get_operation_id_number(value):
    """
    Validates then extracts from PAPI operation IDs just the final number.
    eg:
        'projects/project_name/operations/01234567891011121314' -> '01234567891011121314'
    """
    m = operation_regex.search(value)
    if m:
        return m.group(1)
    else:
        msg = f'Unexpected operation ID {value}. Expected something like {operation_regex.pattern}'
        raise Exception(msg)


def determine_papi_version_from_operation_id(value):
    return Exception("Not Implemented")
