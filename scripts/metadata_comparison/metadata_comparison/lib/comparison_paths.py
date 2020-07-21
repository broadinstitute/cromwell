from google.cloud import storage
import google.auth
import logging
from metadata_comparison.lib import argument_regex

from pathlib import Path, PosixPath
from typing import AnyStr, Union
from abc import ABC, abstractmethod


class ComparisonPath(ABC):
    """
    Abstract Base Class for Local and GCS paths sharing an interface for the purpose of PAPI metadata comparison.
    There's nothing particularly "Comparison" about these paths, I just couldn't think of a better name.
    """
    @staticmethod
    def create(path: Union[AnyStr, Path]):
        if isinstance(path, PosixPath):
            return LocalPath(path)
        elif path.startswith('gs://'):
            bucket, obj = argument_regex.gcs_path_regex_validator(path)
            return GcsPath(bucket, obj)
        return LocalPath(path)

    @staticmethod
    def is_valid_path_string(path: AnyStr) -> bool:
        # ick
        return GcsPath.is_valid_path_string(path) or LocalPath.is_valid_path_string(path)

    @abstractmethod
    def read_text(self, encoding: AnyStr = 'utf_8') -> AnyStr: pass

    # `/` operator, used to implement pathlib.Path style `<existing path> / <new path element>` syntax.
    @abstractmethod
    def __truediv__(self, other): pass

    @abstractmethod
    def exists(self) -> bool: pass

    @abstractmethod
    def mkdir_p(self) -> None: pass

    @abstractmethod
    def write_text(self, content: AnyStr, encoding: AnyStr = 'utf_8') -> None: pass

    @abstractmethod
    def description(self) -> AnyStr: pass


class GcsPath(ComparisonPath):
    def __init__(self, bucket: AnyStr, obj: AnyStr, storage_bucket: storage.Bucket = None):
        self._bucket = bucket
        self._object = obj
        self._storage_blob = None
        if storage_bucket is None:
            credentials, project_id = google.auth.default()
            logging.info(f'Creating storage client for bucket {bucket}')
            client = storage.Client(credentials=credentials)
            self._storage_bucket = client.bucket(bucket)
        else:
            self._storage_bucket = storage_bucket

    def __storage_blob(self) -> storage.Blob:
        if self._storage_blob is None:
            logging.info(f'Creating storage blob for {self}')
            self._storage_blob = self._storage_bucket.blob(self._object)
        return self._storage_blob

    def read_text(self, encoding: AnyStr = 'utf_8') -> AnyStr:
        return self.__storage_blob().download_as_string()

    def __truediv__(self, other) -> ComparisonPath:
        return GcsPath(bucket=self._bucket,
                       obj=f'{self._object}/{other}',
                       storage_bucket=self._storage_bucket)

    def exists(self) -> bool:
        return self.__storage_blob().exists()

    def mkdir_p(self) -> None:
        # Nothing to do here, "directory structure" is implicitly "mkdir -p"'d in GCS.
        pass

    def write_text(self, content: AnyStr, encoding: AnyStr = 'utf_8') -> None:
        self.__storage_blob().upload_from_string(content)

    @staticmethod
    def is_valid_path_string(path: AnyStr) -> bool:
        if path.startswith('gs://'):
            return argument_regex.gcs_path_regex_validator(path)
        return False

    def __str__(self) -> AnyStr:
        return f'gs://{self._bucket}/{self._object}'

    def description(self) -> AnyStr:
        return 'GCS'


class LocalPath(ComparisonPath):
    def __init__(self, local_spec: Union[AnyStr, Path]):
        self.path = Path(local_spec)

    def read_text(self, encoding: AnyStr = 'utf_8') -> AnyStr:
        return self.path.read_text(encoding)

    def __truediv__(self, other) -> ComparisonPath:
        return LocalPath(self.path / other)

    def exists(self) -> bool:
        return self.path.exists()

    def mkdir_p(self) -> None:
        self.path.mkdir(parents=True, exist_ok=True)

    def write_text(self, content: AnyStr, encoding: AnyStr = 'utf_8') -> None:
        self.path.write_text(content, encoding)

    @staticmethod
    def is_valid_path_string(path: AnyStr) -> bool:
        return True

    def __str__(self) -> AnyStr:
        return str(self.path)

    def description(self) -> AnyStr:
        return 'Local filesystem'


def validate_path(p: AnyStr) -> AnyStr:
    if ComparisonPath.is_valid_path_string(p):
        return p
    raise ValueError(f'{p} is not a valid path whatsoever')
