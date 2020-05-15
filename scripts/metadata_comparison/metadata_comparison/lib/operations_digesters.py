
from abc import ABC, abstractmethod
import dateutil.parser
from metadata_comparison.lib.operation_ids import JsonObject, operation_id_to_api_version, \
    PAPI_V1_API_VERSION, PAPI_V2_ALPHA1_API_VERSION, PAPI_V2_BETA_API_VERSION
from datetime import datetime
import re

from typing import AnyStr, Iterator


class OperationDigester(ABC):
    """
    Abstract Base Class for PAPI operation subclasses sharing an interface for the purpose of treating digesters
    uniformly regardless of PAPI version.
    """

    def __init__(self, operation_json: JsonObject):
        self.__json = operation_json

    def __metadata(self) -> JsonObject:
        return self.__json.get('metadata')

    def __events(self) -> JsonObject:
        return self.__metadata()['events']

    def start_time(self) -> AnyStr:
        return self.__metadata().get('createTime')

    def end_time(self) -> AnyStr:
        return self.__metadata().get('endTime')

    def total_time_seconds(self) -> float:
        return (dateutil.parser.parse(self.end_time()) - dateutil.parser.parse(self.start_time())).total_seconds()

    @staticmethod
    def create(operation_json: JsonObject):
        operation_id = operation_json.get('name')
        version = operation_id_to_api_version(operation_id)
        if version == PAPI_V1_API_VERSION:
            return PapiV1OperationDigester(operation_json)
        elif version == PAPI_V2_ALPHA1_API_VERSION:
            return PapiV2AlphaOperationDigester(operation_json)
        elif version == PAPI_V2_BETA_API_VERSION:
            return PapiV2BetaOperationDigester(operation_json)
        else:
            raise ValueError(f"Unrecognized format for PAPI operation ID {operation_id}")

    @abstractmethod
    def docker_image_pull_seconds(self) -> float: pass

    def event_with_description(self, description: AnyStr) -> JsonObject:
        def has_description(event: JsonObject) -> bool:
            return event.get('description') == description

        for unique in filter(has_description, self.__metadata().get('events')):
            return unique

    def event_with_description_like(self, description: AnyStr) -> Iterator[JsonObject]:
        regex = re.compile(description)

        def has_description_like(event: JsonObject) -> bool:
            return regex.match(event.get('description')) is not None

        return filter(has_description_like, self.__metadata().get('events'))


class PapiV1OperationDigester(OperationDigester):
    def __init__(self, operation_json: JsonObject):
        super(PapiV1OperationDigester, self).__init__(operation_json)

    def docker_image_pull_seconds(self) -> float:
        descriptions = ['localizing-files', 'pulling-image']
        end, start = [dateutil.parser.parse(self.event_with_description(d).get('startTime')) for d in descriptions]
        return (end - start).total_seconds()


class PapiV2OperationDigester(OperationDigester, ABC):
    def __init__(self, operation_json: JsonObject):
        super(PapiV2OperationDigester, self).__init__(operation_json)

    def docker_image_pull_seconds(self) -> float:
        description = "^(Started|Stopped) pulling .*"
        pull_events = [dateutil.parser.parse(d.get('timestamp')) for d in self.event_with_description_like(description)]
        pull_events.sort()
        return (pull_events[-1] - pull_events[0]).total_seconds()


class PapiV2AlphaOperationDigester(PapiV2OperationDigester):
    def __init__(self, operation_json: JsonObject):
        super(PapiV2AlphaOperationDigester, self).__init__(operation_json)


class PapiV2BetaOperationDigester(PapiV2OperationDigester):
    def __init__(self, operation_json: JsonObject):
        super(PapiV2BetaOperationDigester, self).__init__(operation_json)
