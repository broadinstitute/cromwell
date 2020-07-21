#!/usr/bin/env python3

import re
from typing import Any, AnyStr, Callable, Dict, List, Sequence, TypeVar, Union

PAPI_V1_OPERATION_REGEX = re.compile('^operations/[^/]*')
PAPI_V2ALPHA1_OPERATION_REGEX = re.compile('^projects/[^/]*/operations/[0-9]*')
PAPI_V2BETA_OPERATION_REGEX = re.compile('^projects/[^/]*/locations/[^/]*/operations/[0-9]*')

PAPI_V1_API_VERSION = 'v1alpha2'
PAPI_V2_ALPHA1_API_VERSION = 'v2alpha1'
PAPI_V2_BETA_API_VERSION = 'v2beta'


def get_operation_id_number(value: str) -> str:
    """
    Validates then extracts from PAPI operation IDs just the final number.
    eg:
        papiv1:       'operations/EMj9o52aLhj78ZLxzunkiHcg0e2BmaAdKg9wcm9kdWN0aW9uUXVldWU -> EMj9o52aLhj78ZLxzunkiHcg0e2BmaAdKg9wcm9kdWN0aW9uUXVldWU'
        papiv2alpha1: 'projects/project_name/operations/01234567891011121314' -> '01234567891011121314'
    """
    return value.split('/')[-1]


def operation_id_to_api_version(value: str) -> str:
    """
    Examines an operation ID and returns the PAPI API version which produced it
    Luckily, this is currently a 1:1 format-to-api mapping so we don't need any other clues to tell the API version.
    """
    if PAPI_V1_OPERATION_REGEX.match(value):
        return PAPI_V1_API_VERSION
    elif PAPI_V2ALPHA1_OPERATION_REGEX.match(value):
        return PAPI_V2_ALPHA1_API_VERSION
    elif PAPI_V2BETA_OPERATION_REGEX.match(value):
        return PAPI_V2_BETA_API_VERSION
    else:
        raise Exception(f'Cannot deduce PAPI api version from unexpected operation ID format \'{value}\'')


# What a JSON Object works out to be.
JsonObject = Dict[str, Union[Union[None, AnyStr, float], Any]]
Accumulator = TypeVar('Accumulator')
OperationId = AnyStr
CallNameSequence = Sequence[AnyStr]

OperationMappingCallFunction = Callable[[Accumulator, OperationId, CallNameSequence, JsonObject], None]


def visit_papi_operations(json_metadata: JsonObject,
                          call_fn: OperationMappingCallFunction,
                          initial_accumulator: Accumulator) -> Accumulator:
    """
    Visits all PAPI operations represented in the Cromwell metadata of `json_metadata`.
    There will be more operations than calls if any of the operations were preempted or were failed and retried.
    For every PAPI operation, the function `call_fn` is invoked with the parameters:

    - Accumulator: an object of the type specified by the caller of this function per `initial_accumulator`.
    - OperationId: The PAPI operation ID of the job as a string.
    - CallNameSequence: The "breadcrumbs" leading to this call (e.g. [grandparent_wf, parent_wf, wf, call])
    - JsonObject: The JSON object representing the individual job being examined.

    The final Accumulator is returned as the result of this function.
    """

    accumulator = initial_accumulator

    def examine_calls(calls: JsonObject, path_so_far: List[AnyStr]) -> None:
        for call_name in calls:
            attempts = calls[call_name]
            for attempt in attempts:
                operation_id = attempt.get('jobId')
                sub_workflow_metadata = attempt.get('subWorkflowMetadata')
                path = build_call_path(call_name, path_so_far, attempt)
                if operation_id:
                    call_fn(accumulator, operation_id, path, attempt)
                if sub_workflow_metadata:
                    examine_calls(sub_workflow_metadata.get('calls', {}), path)

    def build_call_path(call_name: str, path_so_far: List[AnyStr], attempt: dict) -> List[AnyStr]:
        call_path = path_so_far.copy()

        # Remove confusing duplication in subworkflow call names.
        # A parent workflow would name a subworkflow call "parent_wf.sub_wf".
        # The subworkflow would name its calls "sub_wf.sub_call".
        # If those call components were simply joined the result would be
        # "parent_wf.sub_wf.sub_wf.sub_call". This logic removes the duplication of "sub_wf",
        # resulting in "parent_wf.sub_wf.sub_call".
        deduplicated_call_name = call_name
        if len(path_so_far) > 0:
            this_call_components = call_name.split('.')
            if len(this_call_components) > 1 and path_so_far[-1].endswith('.' + this_call_components[0]):
                deduplicated_call_name = '.'.join(this_call_components[1:])

        call_path.append(deduplicated_call_name)
        shard_index = attempt.get('shardIndex', -1)
        if shard_index != -1:
            call_path.append(f"shard_{shard_index:04d}")

        return call_path

    examine_calls(calls=json_metadata.get('calls', {}), path_so_far=[])

    return accumulator
