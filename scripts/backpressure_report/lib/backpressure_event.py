
from datetime import datetime
from typing import AnyStr


class BackpressureEvent:
    def __init__(self, pod: str, start: datetime, end: datetime):
        self.pod = pod
        self.start = start
        self.end = end

    def duration(self):
        return (self.end - self.start).seconds

    def __str__(self) -> AnyStr:
        return f"BackpressureEvent(pod = {self.pod},start={str(self.start)},duration={(self.duration())}s)"
