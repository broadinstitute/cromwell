# from google.cloud import storage
# import google.auth
# import logging as log
import metadata_comparison.lib.argument_regex as argument_regex

from pathlib import Path
from typing import AnyStr, Tuple, Union
from abc import ABC, abstractmethod

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
            raise ValueError(f'Unsupported as DigesterPath: {path}')

    @abstractmethod
    def read_text(self, encoding: AnyStr = 'utf_8') -> AnyStr:
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
    def write_text(self, content: AnyStr, encoding: AnyStr = 'utf_8') -> None:
        pass

    @staticmethod
    def is_valid_path_string(path: AnyStr) -> bool:
        # ick
        return GcsPath.is_valid_path_string(path) or LocalPath.is_valid_path_string(path)


class GcsPath(DigesterPath):
    def __init__(self, bucket: GcsBucket, obj: GcsObject):
        self._bucket = bucket
        self._object = obj

    def read_text(self, encoding: AnyStr = 'utf_8') -> AnyStr:
        raise ValueError("implement me")

    def __truediv__(self, other):
        return GcsPath(self._bucket, '/'.join((Path(self._object) / other).parts))

    def exists(self) -> bool:
        raise ValueError("implement me")

    def mkdir(self) -> None:
        # Nothing to do here, "directory structure" is implicitly "mkdir -p"'d in GCS.
        pass

    def write_text(self, content: AnyStr, encoding: AnyStr = 'utf_8') -> None:
        raise ValueError("implement me")

    @staticmethod
    def is_valid_path_string(path: AnyStr) -> bool:
        if path.startswith('gs://'):
            return argument_regex.gcs_path_regex_validator(path)
        return False


class LocalPath(DigesterPath):
    def __init__(self, local_spec: Union[LocalSpec, Path]):
        self.path = Path(local_spec)

    def read_text(self, encoding: AnyStr = 'utf_8') -> AnyStr:
        return self.path.read_text(encoding)

    def __truediv__(self, other):
        return LocalPath(self.path / other)

    def exists(self) -> bool:
        return self.path.exists()

    def mkdir(self) -> None:
        self.path.mkdir(parents=True, exist_ok=True)

    def write_text(self, content: AnyStr, encoding: AnyStr = 'utf_8') -> None:
        self.path.write_text(content, encoding)

    @staticmethod
    def is_valid_path_string(path: AnyStr) -> bool:
        return True
