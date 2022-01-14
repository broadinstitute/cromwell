
from dateutil import parser

class BackpressureWindow:
    def __init__(self, pod: str, start: str, end: str):
        self.pod = pod
        self.start = parser.isoparse(start)
        self.end = parser.isoparse(end)
