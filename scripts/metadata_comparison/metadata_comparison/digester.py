import argparse
import json
from google.cloud import storage
import google.auth
# from googleapiclient.discovery import build as google_client_build
# import logging as log
import metadata_comparison.lib.argument_regex as argument_regex
import metadata_comparison.lib.operation_ids as operation_ids
from metadata_comparison.lib.operation_ids import Accumulator, CallNameSequence, JsonObject, OperationId

# from datetime import timedelta, datetime, timezone
import dateutil.parser
import os
from pathlib import Path
from typing import Any, AnyStr, Dict, List, Tuple, Union
from abc import ABC, abstractmethod

Version = "0.0.1"
verbose = False


GcsBucket = AnyStr
GcsObject = AnyStr
GcsSpec = Tuple[GcsBucket, GcsObject]
LocalSpec = Union[AnyStr, Path]


class DigesterPath(ABC):
    @staticmethod
    def create(path: Union[GcsSpec, LocalSpec]):
        if isinstance(path, Tuple):
            return GcsPath(path[0], path[1])
        elif isinstance(path, bytes) or isinstance(path, str) or isinstance(path, Path):
            return LocalPath(path)
        else:
            raise ValueError(f'Unrecognized path {path}')

    @abstractmethod
    def read_bytes(self) -> AnyStr:
        pass

    @abstractmethod
    def __truediv__(self, other) -> AnyStr:
        pass

    @abstractmethod
    def exists(self) -> bool:
        pass

    @abstractmethod
    def mkdir(self) -> None:
        pass

    @abstractmethod
    def write_bytes(self, content: AnyStr) -> None:
        pass


class GcsPath(DigesterPath):
    def __init__(self, bucket: GcsBucket, obj: GcsObject):
        self._bucket = bucket
        self._object = obj

    def read_bytes(self) -> AnyStr:
        raise ValueError("implement me")

    def __truediv__(self, other):
        return GcsPath(self._bucket, '/'.join((Path(self._object) / other).parts))

    def exists(self) -> bool:
        raise ValueError("implement me")

    def mkdir(self) -> None:
        # Nothing to do here, "directory structure" is implicitly "mkdir -p"'d in GCS.
        pass

    def write_bytes(self, content: AnyStr) -> None:
        raise ValueError("implement me")


class LocalPath(DigesterPath):
    def __init__(self, local_spec: Union[LocalSpec, Path]):
        self.path = Path(local_spec)

    def read_bytes(self) -> AnyStr:
        return self.path.read_bytes()

    def __truediv__(self, other):
        return LocalPath(self.path / other)

    def exists(self) -> bool:
        return self.path.exists()

    def mkdir(self) -> None:
        self.path.mkdir(parents=True, exist_ok=True)

    def write_bytes(self, content: AnyStr) -> None:
        self.path.write_bytes()


def main() -> None:
    args = parse_args()

    for path in args.local_paths:
        parent_path = DigesterPath.create(path)

        workflow_path = parent_path / 'workflow.json'
        operations_dir_path = parent_path / 'operations'

        digest_parent = parent_path / 'digests' / Version
        digest_path = digest_parent / 'digest.json'

        if not digest_path.exists() or args.force:
            digest_parent.mkdir()
            digest_json = digest(workflow_path, operations_dir_path)
            digest_string = json.dumps(digest_json, sort_keys=True, indent=4)
            digest_path.write(digest_string)
        else:
            raise ValueError(f'digest file already exists at {digest_path} and --force not specified')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Digest workflow metadata and job operation details, reading from and reuploading to GCS.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='whether to log verbosely (default False)')
    parser.add_argument('-f', '--force', action='store_true',
                        help='whether to overwrite existing digests (default False)')
    parser.add_argument('local_paths', metavar="PATH", nargs='+', type=argument_regex.validate_gcs_if_gcs,
                        help="Location at which to find metadata (local or GCS)")

    return parser.parse_args()


CallName = AnyStr


def digest(workflow_path: DigesterPath, operations_path: DigesterPath) -> JsonObject:
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

            succeeded_operations[string_path] = {
                "attempt": attempt.get('attempt'),
                "shardIndex": attempt.get('shardIndex'),
                "operationId": operation_id,
                "cromwellStart": cromwell_start,
                "cromwellEnd": cromwell_end,
                "cromwellTotalTimeSeconds": cromwell_total_time_seconds
            }

    data = workflow_path.read_bytes()
    metadata = json.loads(data)

    shards = operation_ids.visit_papi_operations(metadata, call_fn, initial_accumulator={})
    return {'version': Version, 'calls': shards, 'workflowId': metadata['id']}


if __name__ == "__main__":
    main()
