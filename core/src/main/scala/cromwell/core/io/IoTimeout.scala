package cromwell.core.io

import java.util.concurrent.TimeoutException

case class IoTimeout(command: IoCommand[_]) extends TimeoutException(s"The I/O operation $command timed out")
