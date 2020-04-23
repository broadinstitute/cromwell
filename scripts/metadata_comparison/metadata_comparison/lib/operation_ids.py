#!/usr/bin/env python3

import re

papi_v1_operation_regex = re.compile('^operations/[^/]*')
papi_v2alpha1_operation_regex = re.compile('^projects/.*/operations/[0-9]*')
papi_v2beta_operation_regex = re.compile('^projects/.*/locations/.*/operations/[0-9]*')

def get_operation_id_number(value):
    """
    Validates then extracts from PAPI operation IDs just the final number.
    eg:
        papiv1:       'operations/EMj9o52aLhj78ZLxzunkiHcg0e2BmaAdKg9wcm9kdWN0aW9uUXVldWU -> EMj9o52aLhj78ZLxzunkiHcg0e2BmaAdKg9wcm9kdWN0aW9uUXVldWU'
        papiv2alpha1: 'projects/project_name/operations/01234567891011121314' -> '01234567891011121314'
    """
    return value.split('/')[-1]


def operation_id_to_api_version(value):
    """
    Examines an operation ID and returns the PAPI API version which produced it
    Luckily, this is currently a 1:1 format-to-api mapping so we don't need any other clues to tell the API version.
    """
    if papi_v1_operation_regex.match(value):
        return 'v1alpha2'
    elif papi_v2alpha1_operation_regex.match(value):
        return 'v2alpha1'
    elif papi_v2beta_operation_regex.match(value):
        return 'v2beta'
    else:
        raise Exception(f'Cannot deduce PAPI api version from ID "{value}"')


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
