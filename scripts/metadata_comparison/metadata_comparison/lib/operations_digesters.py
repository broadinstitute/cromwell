
from abc import ABC, abstractmethod
import dateutil.parser
from metadata_comparison.lib.operation_ids import JsonObject, operation_id_to_api_version, \
    PAPI_V1_API_VERSION, PAPI_V2_ALPHA1_API_VERSION, PAPI_V2_BETA_API_VERSION
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

    def create_time(self) -> AnyStr:
        return self.__metadata().get('createTime')

    def start_time(self) -> AnyStr:
        return self.__metadata().get('startTime')

    def end_time(self) -> AnyStr:
        return self.__metadata().get('endTime')

    def total_time_seconds(self) -> float:
        return (dateutil.parser.parse(self.end_time()) - dateutil.parser.parse(self.create_time())).total_seconds()

    def metadata(self):
        return self.__metadata()

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
    def docker_image_pull_time_seconds(self) -> float: pass

    @abstractmethod
    def localization_time_seconds(self) -> float: pass

    @abstractmethod
    def user_command_time_seconds(self) -> float: pass

    @abstractmethod
    def delocalization_time_seconds(self) -> float: pass

    @abstractmethod
    def startup_time_seconds(self) -> float: pass

    @abstractmethod
    def machine_type(self) -> AnyStr: pass

    def other_time_seconds(self) -> float:
        end, create = [dateutil.parser.parse(t) for t in [self.end_time(), self.create_time()]]
        total_time = (end - create).total_seconds()

        accounted_for_time = \
            self.startup_time_seconds() + \
            self.docker_image_pull_time_seconds() + \
            self.localization_time_seconds() + \
            self.user_command_time_seconds() + \
            self.delocalization_time_seconds()

        return max(total_time - accounted_for_time, 0)

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

    @staticmethod
    def __total_seconds_between_timestamps(start: AnyStr, end: AnyStr) -> float:
        _start, _end = [dateutil.parser.parse(t) for t in [start, end]]
        return (_end - _start).total_seconds()

    @staticmethod
    def __total_seconds_between_events(start: JsonObject, end: JsonObject) -> float:
        _start, _end = [e.get('startTime') for e in [start, end]]
        return PapiV1OperationDigester.__total_seconds_between_timestamps(_start, _end)

    def startup_time_seconds(self) -> float:
        # Look at `pulling_image` as that is the next lifecycle phase after startup.
        create, pulling_image = self.create_time(), self.event_with_description('pulling-image').get('startTime')
        return self.__total_seconds_between_timestamps(create, pulling_image)

    def docker_image_pull_time_seconds(self) -> float:
        start, end = [self.event_with_description(d) for d in ['pulling-image', 'localizing-files']]
        return self.__total_seconds_between_events(start, end)

    def localization_time_seconds(self) -> float:
        start, end = [self.event_with_description(d) for d in ['localizing-files', 'running-docker']]
        return self.__total_seconds_between_events(start, end)

    def user_command_time_seconds(self) -> float:
        start, end = [self.event_with_description(d) for d in ['running-docker', 'delocalizing-files']]
        return self.__total_seconds_between_events(start, end)

    def delocalization_time_seconds(self) -> float:
        start, end = [self.event_with_description(d) for d in ['delocalizing-files', 'ok']]
        return self.__total_seconds_between_events(start, end)

    def machine_type(self) -> AnyStr:
        machine_type_with_zone_prefix = self.metadata().get('runtimeMetadata').get('computeEngine').get('machineType')
        return machine_type_with_zone_prefix.split('/')[-1]


class PapiV2OperationDigester(OperationDigester, ABC):
    def __init__(self, operation_json: JsonObject):
        super(PapiV2OperationDigester, self).__init__(operation_json)

    def startup_time_seconds(self) -> float:
        create = dateutil.parser.parse(self.create_time())
        docker_description = "^Started pulling .*"
        docker_events = [dateutil.parser.parse(d.get('timestamp')) for d in
                         self.event_with_description_like(docker_description)]
        docker_events.sort()
        return (docker_events[0] - create).total_seconds()

    def docker_image_pull_time_seconds(self) -> float:
        description = "^(Started|Stopped) pulling .*"
        events = [dateutil.parser.parse(d.get('timestamp')) for d in self.event_with_description_like(description)]
        events.sort()
        return (events[-1] - events[0]).total_seconds()

    def localization_time_seconds(self) -> float:
        description = "^.* (Starting|Done)\\\\\\\\ localization.\"$"
        events = [dateutil.parser.parse(d.get('timestamp')) for d in self.event_with_description_like(description)]
        events.sort()
        return (events[-1] - events[0]).total_seconds()

    def user_command_time_seconds(self) -> float:
        started_running_description = "^Started running \"/cromwell_root/script\"$"
        started_running_events = [dateutil.parser.parse(d.get('timestamp')) for d in
                                  self.event_with_description_like(started_running_description)]

        stopped_running_description = "^Stopped running \"/cromwell_root/script\"$"
        stopped_running_events = [dateutil.parser.parse(d.get('timestamp')) for d in
                                  self.event_with_description_like(stopped_running_description)]

        events = started_running_events + stopped_running_events
        events.sort()
        return (events[-1] - events[0]).total_seconds()

    def delocalization_time_seconds(self) -> float:
        description = "^.* (Starting|Done)\\\\\\\\ delocalization.\"$"
        events = [dateutil.parser.parse(d.get('timestamp')) for d in self.event_with_description_like(description)]
        events.sort()
        return (events[-1] - events[0]).total_seconds()

    def machine_type(self) -> AnyStr:
        event = next(self.event_with_description_like('^Worker .* assigned in .*'))
        return event.get('details').get('machineType')


class PapiV2AlphaOperationDigester(PapiV2OperationDigester):
    def __init__(self, operation_json: JsonObject):
        super(PapiV2AlphaOperationDigester, self).__init__(operation_json)


class PapiV2BetaOperationDigester(PapiV2OperationDigester):
    def __init__(self, operation_json: JsonObject):
        super(PapiV2BetaOperationDigester, self).__init__(operation_json)
