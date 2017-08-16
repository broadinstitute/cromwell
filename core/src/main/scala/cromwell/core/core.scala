package cromwell.core

import cromwell.core.path.Path
import lenthall.exception.ThrowableAggregation
import wdl4s.wdl.values.WdlValue


case class CallContext(root: Path, stdout: String, stderr: String)
case class JobOutput(wdlValue: WdlValue)
/**  Marker trait for Cromwell exceptions that are to be treated as fatal (non-retryable) */
trait CromwellFatalExceptionMarker { this: Throwable => }

object CromwellFatalException {
  // Don't wrap if it's already a fatal exception
  def apply(throwable: Throwable) = throwable match {
    case e: CromwellFatalExceptionMarker => e
    case e => new CromwellFatalException(e)
  }
  def unapply(e: CromwellFatalException): Option[Throwable] = Option(e.exception)
}

class CromwellFatalException(val exception: Throwable) extends Exception(exception) with CromwellFatalExceptionMarker

case class CromwellAggregatedException(throwables: Seq[Throwable], exceptionContext: String = "")
  extends Exception with ThrowableAggregation
