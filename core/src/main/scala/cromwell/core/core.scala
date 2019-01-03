package cromwell.core

import common.exception.ThrowableAggregation
import cromwell.core.path.Path

import scala.concurrent.duration.FiniteDuration


case class StandardPaths(output: Path, error: Path)

case class CallContext(root: Path, standardPaths: StandardPaths, isDocker: Boolean)

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

case class CacheConfig(concurrency: Int, size: Long, ttl: FiniteDuration)
