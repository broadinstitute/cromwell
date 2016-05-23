package cromwell

import lenthall.exception.ThrowableAggregation
import java.nio.file.Path

import wdl4s.values.{SymbolHash, WdlValue}

import scalaz._

package object core {
  // root can be a Path instead of a String in PBE. stdout / err too but it doesn't really bring values since they're just stringified to WdlFiles
  class OldWorkflowContext(val root: String)
  class OldCallContext(override val root: String, val stdout: String, val stderr: String) extends OldWorkflowContext(root)

  case class CallContext(root: Path, stdout: String, stderr: String)

  type ErrorOr[+A] = ValidationNel[String, A]
  type LocallyQualifiedName = String
  case class CallOutput(wdlValue: WdlValue, hash: Option[SymbolHash])
  type JobOutputs = Map[LocallyQualifiedName, CallOutput]
  type EvaluatedRuntimeAttributes = Map[String, WdlValue]

  class CromwellFatalException(exception: Throwable) extends Exception(exception)
  case class CromwellAggregatedException(throwables: Seq[Throwable], exceptionContext: String = "") extends ThrowableAggregation
}
