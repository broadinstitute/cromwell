#!/usr/bin/env python3
#
# comparer.py
#
# Purpose: Compare performance metadata JSON files produced by Digester and produce result in CSV format
#
# Usage: python3 -m metadata_comparison.comparer [-h] [-v] [--force]
#                --name1 NAME_FOR_DIGEST_1 --name2 NAME_FOR_DIGEST_2
#                --digest1 PATH_TO_DIGEST_1 --digest2 PATH_TO_DIGEST_2 --output-path OUTPUT_PATH
#                [--call-prefix-to-remove [CALL_PREFIX_TO_REMOVE [CALL_PREFIX_TO_REMOVE ...]]]
#
# For ExomeGermlineSingleSample workflows the local call names are globally unique so all
# FQN prefixes can be removed for ease of interpretation. An invocation to compare PAPI v1
# to PAPI v2 might look like:
#
# python3 -m metadata_comparison.comparer --name1 PAPIv1 --name2 PAPIv2 \
#         --digest1 papiv1.json --digest2 papiv2.json --output-path comparison.csv \
#         --call-prefix-to-remove ExomeGermlineSingleSample.AggregatedBamQC. \
#         ExomeGermlineSingleSample.BamToCram. ExomeGermlineSingleSample.BamToGvcf.VariantCalling. \
#         ExomeGermlineSingleSample.UnmappedBamToAlignedBam. ExomeGermlineSingleSample.
#
# Python Prereqs (at least, the ones which I needed to manually install... YMMV):
#
#   * pip3 install --upgrade google-api-python-client
#   * pip3 install --upgrade google-cloud-storage
#
# Remember to login to create application default credentials before use:
#   % gcloud auth application-default login

import argparse
import json
import logging
from metadata_comparison.lib.digester_keys import *
from metadata_comparison.lib.comparison_paths import ComparisonPath, validate_path
from metadata_comparison.lib.logging import set_log_verbosity, quieten_chatty_imports
from metadata_comparison.lib.operation_ids import JsonObject
from typing import AnyStr, List

logger = logging.getLogger('metadata_comparison.comparer')


class DigesterKey:
    def __init__(self, json_key: AnyStr, display_text: AnyStr):
        self.json_key = json_key
        self.display_text = display_text


DigesterKeys = [
    DigesterKey(json_key=PapiTotalTimeSeconds, display_text='Total PAPI time'),
    DigesterKey(json_key=StartupTimeSeconds, display_text='Startup'),
    DigesterKey(json_key=DockerImagePullTimeSeconds, display_text='Docker Pull'),
    DigesterKey(json_key=LocalizationTimeSeconds, display_text='Localization'),
    DigesterKey(json_key=UserCommandTimeSeconds, display_text='User command'),
    DigesterKey(json_key=DelocalizationTimeSeconds, display_text='Delocalization'),
    DigesterKey(json_key=OtherTimeSeconds, display_text='Other time'),
    DigesterKey(json_key=MachineType, display_text='Machine type')
]


def digester_key_by_json_key(json_key: AnyStr) -> DigesterKey:
    return next(dk for dk in DigesterKeys if dk.json_key == json_key)


MachineTypesCostPerHour = {
    'n1-highcpu-16': 0.5672,    # https://cloud.google.com/compute/vm-instance-pricing#n1_highcpu_machine_types
    'n1-highmem-2': 0.1184,     # https://cloud.google.com/compute/vm-instance-pricing#n1_highmem_machine_types
    'n1-standard-1': 0.0475,    # https://cloud.google.com/compute/vm-instance-pricing#n1_standard_machine_types
    'n1-standard-2': 0.0950,    # https://cloud.google.com/compute/vm-instance-pricing#n1_standard_machine_types
    'g1-small': 0.0257,         # https://cloud.google.com/compute/all-pricing#n1_sharedcore_machine_types
    'e2-standard-2': 0.067006,  # https://cloud.google.com/compute/vm-instance-pricing#e2_standard_machine-types
    'n2-standard-2': 0.0971     # https://cloud.google.com/compute/vm-instance-pricing#n2_standard_machine_types
}

# A machine type-weighted dictionary used for getting more accurate cost estimates.
MachineTypeCostMultiplier = {}
for key in MachineTypesCostPerHour.keys():
    # Normalize to g1-small, currently the least expensive machine type used.
    MachineTypeCostMultiplier[key] = MachineTypesCostPerHour[key] / MachineTypesCostPerHour['g1-small']


class CallKey:
    def __init__(self, full, call_prefixes_to_remove: List[AnyStr]):
        self.full = full
        # Call keys without repetitive prefixes are useful for building comparer reports.
        # Assign `without_prefix` to the full prefix in case none of the `call_prefixes_to_remove` match.
        self.without_prefix = full
        for prefix in call_prefixes_to_remove:
            if self.full.startswith(prefix):
                self.without_prefix = self.full[len(prefix):]
                break


def format_seconds(total_seconds: float) -> AnyStr:
    """
    Format the specified number of seconds as HH:MM:SS.
    """
    hours, remainder = divmod(int(total_seconds), 3600)
    minutes, seconds = divmod(remainder, 60)
    return f'{hours}:{minutes:02}:{seconds:02}'


def compare_jsons(json_1: JsonObject, json_2: JsonObject,
                  name_1: AnyStr, name_2: AnyStr,
                  call_prefixes_to_remove: List[AnyStr],
                  force: bool = False) -> List[List[AnyStr]]:
    """
    Produce a CSV representing the comparison of the specified JSONs.
    """
    error_checks(name_1, name_2, json_1, json_2, force)

    call_keys = [CallKey(k, call_prefixes_to_remove) for k in json_1.get('calls').keys()]

    call_keys_sorted_without_prefix = sorted(call_keys, key=lambda k: k.without_prefix)

    # Columns to produce in triplets: name_1, name_2, percent increase from name_1 to name_2.
    digester_key_names = [PapiTotalTimeSeconds, StartupTimeSeconds, DockerImagePullTimeSeconds, LocalizationTimeSeconds,
                          UserCommandTimeSeconds, DelocalizationTimeSeconds]

    def build_header_rows():
        """ Builds the header rows of the spreadsheet. """
        top_header_row = ['job', 'Machine type']
        percent_contrib_row = ['% contribution to total run time', '']
        total_row = ['Total', '']
        machine_type_weighted_row = ['Machine-type weighted total (cost-ish)', '']

        def sum_call_times(j: JsonObject, key: AnyStr):
            """ Sums raw call times for a column. """
            times = [j.get('calls').get(ck.full).get(key) for ck in call_keys]
            return sum(times)

        def sum_call_times_weighted(j: JsonObject, key: AnyStr):
            """
            Sums call times for a column weighted by the machine type used as a more accurate estimation
            of cost.
            """
            times = [j.get('calls').get(ck.full).get(key) *
                     MachineTypeCostMultiplier[j.get('calls').get(ck.full).get('machineType')] for ck in call_keys]
            return sum(times)

        # Treat the total time specially to be able to report percentages of time for individual lifecycle phases.
        total_total_time_1, total_total_time_2 = [
            sum_call_times(j, PapiTotalTimeSeconds) for j in [json_1, json_2]]

        for digester_key_name in digester_key_names:
            digester_key = digester_key_by_json_key(digester_key_name)
            top_header_row.append(f'{name_1} {digester_key.display_text}')
            top_header_row.append(f'{name_2} {digester_key.display_text}')
            top_header_row.append('% increase')

            # Total times within a column, not necessarily total time of the job.
            total_time_1, total_time_2 = [
                sum_call_times(j, digester_key_name) for j in [json_1, json_2]]

            percent_contrib_row.append(f'{(total_time_1 / total_total_time_1) * 100:.2f}%')
            percent_contrib_row.append(f'{(total_time_2 / total_total_time_2) * 100:.2f}%')
            percent_contrib_row.append('')

            total_row.append(format_seconds(total_time_1))
            total_row.append(format_seconds(total_time_2))
            total_percent_increase = ((total_time_2 - total_time_1) * 100) / total_time_1
            total_row.append(f'{total_percent_increase :.2f}%')

            weighted_total_time_1, weighted_total_time_2 = [
                sum_call_times_weighted(j, digester_key_name) for j in [json_1, json_2]]

            machine_type_weighted_row.append(format_seconds(weighted_total_time_1))
            machine_type_weighted_row.append(format_seconds(weighted_total_time_2))
            weighted_percent = ((weighted_total_time_2 - weighted_total_time_1) * 100) / weighted_total_time_1
            machine_type_weighted_row.append(f'{weighted_percent:.2f}%')

        return [top_header_row, percent_contrib_row, total_row, machine_type_weighted_row, []]

    rows = []

    # Produce data for individual lifecycle phases of individual calls.
    for call_key in call_keys_sorted_without_prefix:
        row = []
        call_1, call_2 = [j.get('calls').get(call_key.full) for j in [json_1, json_2]]
        row.append(call_key.without_prefix)
        row.append(call_1.get(MachineType))

        for digester_key_name in digester_key_names:
            time_1 = call_1.get(digester_key_name)
            time_2 = call_2.get(digester_key_name)
            row.append(format_seconds(time_1))
            row.append(format_seconds(time_2))
            if time_1:
                row.append(f'{((time_2 - time_1)  * 100) / time_1:.2f}%')
            else:
                row.append('---')

        rows.append(row)

    return build_header_rows() + rows


def error_checks(name_1: AnyStr, name_2: AnyStr, json_1: JsonObject, json_2: JsonObject, force: bool = False):
    version_1, version_2 = [j.get('version') for j in [json_1, json_2]]

    if version_1 != version_2:
        msg = f"Inconsistent digest versions: First JSON digest is {version_1} but second is {version_2}"
        raise ValueError(msg)

    call_keys_1, call_keys_2 = [list(j.get('calls').keys()) for j in [json_1, json_2]]
    call_keys_1.sort()
    call_keys_2.sort()

    if call_keys_1 != call_keys_2:
        in_1 = ', '.join(set(call_keys_1) - set(call_keys_2))
        in_2 = ', '.join(set(call_keys_2) - set(call_keys_1))
        msg_1 = None
        msg_2 = None
        if in_1:
            msg_1 = f"In {name_1} but not {name_2}: {in_1}."
        if in_2:
            msg_2 = f"In {name_1} but not {name_2}: {in_2}."

        raise ValueError('The specified digest files do not have the same call keys. These digests cannot be ' +
                         f'compared and probably are not from the same workflow and sample. {msg_1} {msg_2}')

    for call_key in call_keys_1:
        call_1 = json_1.get('calls').get(call_key)
        call_2 = json_2.get('calls').get(call_key)
        for call in [call_1, call_2]:
            for digester_key in DigesterKeys:
                json_key = digester_key.json_key
                if json_key not in call:
                    if call == call_1:
                        nth = "first"
                    else:
                        nth = "second"
                    raise ValueError(
                        f"In {nth} digest JSON: call '{call_key}' missing required key '{json_key}'")

    discrepancies = []

    for k in call_keys_1:
        machine_type_1 = json_1.get('calls').get(k).get('machineType')
        machine_type_2 = json_2.get('calls').get(k).get('machineType')
        if machine_type_1 != machine_type_2:
            discrepancies.append({
                'key': k,
                'machine_type_1': machine_type_1,
                'machine_type_2': machine_type_2
            })

    if discrepancies:
        string = ', '.join(f'{e["key"]}: {e["machine_type_1"]} vs {e["machine_type_2"]}' for e in discrepancies)
        message = 'The specified digest files unexpectedly contain corresponding jobs that ran with ' + \
                  'different machine types: ' + string
        if force:
            logger.warning(message)
        else:
            raise ValueError(message + '. Specify the --force argument to force comparison anyway.')


def json_from_path_string(path_string: AnyStr) -> JsonObject:
    path = ComparisonPath.create(path_string)
    return json.loads(path.read_text())


def csv_string_from_data(data: List[List[AnyStr]]) -> AnyStr:
    rows = [','.join(row) for row in data]
    return '\n'.join(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare performance digest JSONs and produce CSV result')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('--force', action='store_true',
                        help='Whether to force the comparer to run through certain warning conditions')
    parser.add_argument('--name1', nargs=1, metavar='NAME_FOR_DIGEST_1', type=str,
                        required=True, help='Name to use for the first digest')
    parser.add_argument('--name2', nargs=1, metavar='NAME_FOR_DIGEST_2', type=str,
                        required=True, help='Name to use for the second digest')
    parser.add_argument('--digest1', nargs=1, metavar='PATH_TO_DIGEST_1', type=validate_path,
                        required=True, help='First digest path to compare.')
    parser.add_argument('--digest2', nargs=1, metavar='PATH_TO_DIGEST_2', type=validate_path,
                        required=True, help='Second digest path to compare.')
    parser.add_argument('--output-path', metavar='OUTPUT_PATH', type=validate_path, nargs=1,
                        required=True, help='Path to output CSV file.')
    parser.add_argument('--call-prefix-to-remove', metavar='CALL_PREFIX_TO_REMOVE', type=str, nargs='*',
                        help='Call prefix to remove if present.')

    args = parser.parse_args()
    set_log_verbosity(args.verbose)
    quieten_chatty_imports()
    logger.info("Starting Comparer operation.")

    _json_1, _json_2 = [json_from_path_string(p[0]) for p in [args.digest1, args.digest2]]

    prefixes = [] if not args.call_prefix_to_remove else args.call_prefix_to_remove

    comparison_data = compare_jsons(_json_1, _json_2, args.name1[0], args.name2[0], prefixes, args.force)
    out_path = args.output_path[0]
    out = ComparisonPath.create(out_path)
    if out.exists() and not args.force:
        raise ValueError(f"Specified output file '{out_path}' already exists and --force not specified.")
    out.write_text(csv_string_from_data(comparison_data))

    logger.info('Comparer operation completed successfully.')
