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


def find_operation_ids_in_metadata(json_metadata):
    """Finds all instances of PAPI operations IDs in a workflow"""
    # Eg given:
    # {
    #   "calls": {
    #     "workflow_name.task_name": [
    #       {
    #         "jobId": "projects/broad-dsde-cromwell-dev/operations/01234567891011121314",
    # ...
    #
    # We want to extract "projects/broad-dsde-cromwell-dev/operations/01234567891011121314"
    papi_operations = []

    def find_operation_ids_in_calls(calls):
        for callname in calls:
            attempts = calls[callname]
            for attempt in attempts:
                operation_id = attempt.get('jobId')
                subWorkflowMetadata = attempt.get('subWorkflowMetadata')
                if operation_id:
                    papi_operations.append(operation_id)
                if subWorkflowMetadata:
                    find_operation_ids_in_calls(subWorkflowMetadata.get('calls', {}))

    find_operation_ids_in_calls(json_metadata.get('calls', {}))

    return papi_operations
