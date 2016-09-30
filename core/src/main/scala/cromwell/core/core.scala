package cromwell.core

import java.nio.file.Path

import lenthall.exception.ThrowableAggregation
import wdl4s.values.WdlValue


case class CallContext(root: Path, stdout: String, stderr: String)
case class JobOutput(wdlValue: WdlValue)
class CromwellFatalException(exception: Throwable) extends Exception(exception)
case class CromwellAggregatedException(throwables: Seq[Throwable], exceptionContext: String = "") extends ThrowableAggregation
