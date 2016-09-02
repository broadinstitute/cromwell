package cromwell

import java.nio.file.Path

import lenthall.exception.ThrowableAggregation
import wdl4s.values.WdlValue

package object core {
  case class CallContext(root: Path, stdout: String, stderr: String)

  type LocallyQualifiedName = String
  type FullyQualifiedName = String
  type WorkflowOutputs = Map[FullyQualifiedName, JobOutput]
  type WorkflowOptionsJson = String
  case class JobOutput(wdlValue: WdlValue)
  type JobOutputs = Map[LocallyQualifiedName, JobOutput]
  type HostInputs = Map[String, WdlValue]
  type EvaluatedRuntimeAttributes = Map[String, WdlValue]

  class CromwellFatalException(exception: Throwable) extends Exception(exception)
  case class CromwellAggregatedException(throwables: Seq[Throwable], exceptionContext: String = "") extends ThrowableAggregation
}
