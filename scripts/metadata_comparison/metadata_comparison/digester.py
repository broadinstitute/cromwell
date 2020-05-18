import argparse
import json
from metadata_comparison.lib import logging, operation_ids
from metadata_comparison.lib.operation_ids import CallNameSequence, JsonObject, OperationId
from metadata_comparison.lib.comparison_paths import ComparisonPath
from metadata_comparison.lib.operations_digesters import OperationDigester

import dateutil.parser
from typing import AnyStr, Dict

Version = "0.0.1"


def main(args: argparse.Namespace) -> None:
    for path in args.paths:
        parent_path = ComparisonPath.create(path)

        workflow_path = parent_path / 'workflow.json'
        operations_dir_path = parent_path / 'operations'

        digest_parent = parent_path / 'digests' / Version
        digest_path = digest_parent / 'digest.json'

        if not digest_path.exists() or args.force:
            digest_parent.mkdir_p()
            digest_json = digest(workflow_path, operations_dir_path)
            digest_string = json.dumps(digest_json, sort_keys=True, indent=4)
            digest_path.write_text(digest_string)
        else:
            raise ValueError(f'digest file already exists at {digest_path} and --force not specified')


def parse_args() -> argparse.Namespace:
    def validate_path(p: AnyStr) -> AnyStr:
        if ComparisonPath.is_valid_path_string(p):
            return p
        raise ValueError(f'{p} is not a valid path whatsoever')

    parser = argparse.ArgumentParser(
        description='Digest workflow metadata and job operation details, reading from and reuploading to GCS.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='whether to log verbosely (default False)')
    parser.add_argument('-f', '--force', action='store_true',
                        help='whether to overwrite existing digests (default False)')
    parser.add_argument('paths', metavar="PATH", nargs='+', type=validate_path,
                        help="Location at which to find metadata (local or GCS)")

    return parser.parse_args()


CallName = AnyStr


def digest(workflow_path: ComparisonPath, operations_path: ComparisonPath) -> JsonObject:
    def call_fn(succeeded_operations: Dict[CallName, JsonObject],
                operation_id: OperationId,
                path: CallNameSequence,
                attempt: JsonObject) -> None:
        backend_status = attempt.get('backendStatus', 'Unknown')
        # This script should only ever be pointed at successful workflow metadata. All jobs that have a backend status
        # other than `Success` must have later been re-run successfully, so any un`Success`ful attempts are ignored.
        # It's possible that a future version of the digester might actually want to look at these jobs since they
        # may have completed some lifecycle events which could be useful in accumulating more performance data.
        if backend_status == 'Success':
            string_path = '.'.join(path)
            cromwell_start = attempt.get('start')
            cromwell_end = attempt.get('end')

            cromwell_total_time_seconds = (dateutil.parser.parse(cromwell_end) -
                                           dateutil.parser.parse(cromwell_start)).total_seconds()

            bare_operation_id = operation_id.split('/')[-1]
            operations_file_path = operations_path / f'{bare_operation_id}.json'
            operations_data = operations_file_path.read_text()
            operations_metadata = json.loads(operations_data)
            operation = OperationDigester.create(operations_metadata)

            papi_total_time_seconds = operation.total_time_seconds()

            cromwell_additional_total_time_seconds = \
                float("%.3f" % (cromwell_total_time_seconds - papi_total_time_seconds))

            succeeded_operations[string_path] = {
                "attempt": attempt.get('attempt'),
                "shardIndex": attempt.get('shardIndex'),
                "operationId": operation_id,
                "cromwellStart": cromwell_start,
                "cromwellEnd": cromwell_end,
                "cromwellTotalTimeSeconds": cromwell_total_time_seconds,
                "papiStart": operation.start_time(),
                "papiEnd": operation.end_time(),
                "papiTotalTimeSeconds": operation.total_time_seconds(),
                "cromwellAdditionalTotalTimeSeconds": cromwell_additional_total_time_seconds,
                "dockerImagePullSeconds": operation.docker_image_pull_seconds()
            }

    data = workflow_path.read_text()
    metadata = json.loads(data)

    shards = operation_ids.visit_papi_operations(metadata, call_fn, initial_accumulator={})
    return {'version': Version, 'calls': shards, 'workflowId': metadata['id']}


if __name__ == "__main__":
    logging.quieten_chatty_imports()
    _args = parse_args()
    logging.set_log_verbosity(_args.verbose)

    main(_args)
