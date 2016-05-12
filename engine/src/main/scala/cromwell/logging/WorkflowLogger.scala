package cromwell.logging

import akka.event.LoggingAdapter
import cromwell.engine.backend.OldStyleWorkflowDescriptor
import org.apache.commons.lang3.exception.ExceptionUtils
import org.slf4j.Logger

/**
 * This sets up two loggers for a given workflow, specified by the
 * WorkflowDescriptor.
 *
 * 1) akkaLogger - the Akka event logger (optional).  The log
 *    messages sent to the Akka logger will propagate to the root
 *    logger (see logback.xml)
 *
 * 2) A logger that writes to a file that only has log messages
 *    for this particular workflow.  This makes it easier to debug
 *    a particular workflow if the messages aren't interleaved with
 *    other workflows.  This logger does NOT propagate to the root
 *    logger
 *
 * 3) otherLoggers - any org.slf4j.Logger instances to also log to.
 *    (defaults to Seq.empty)
 *
 * Each log method called (e.g. warn(), error()) will log to all the
 * above locations
 */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class WorkflowLogger(caller: String,
                          descriptor: OldStyleWorkflowDescriptor,
                          akkaLogger: Option[LoggingAdapter] = None,
                          otherLoggers: Seq[Logger] = Seq.empty[Logger],
                          callTag: Option[String] = None) {
  lazy val workflowFileLogger = descriptor.workflowLogger

  lazy val tag: String = {
    val subtag = callTag.map(c => s":$c").getOrElse("")
    s"$caller [UUID(${descriptor.id.shortString})$subtag]"
  }

  private def format(msg: String): String = s"$tag: $msg"

  private def format(msg: String, throwable: Throwable): String = {
    format(msg) + ":\n" + ExceptionUtils.getStackTrace(throwable)
  }

  def warn(msg: String): Unit = {
    workflowFileLogger.warn(format(msg))
    akkaLogger.foreach(_.warning(format(msg)))
    otherLoggers.foreach(_.warn(format(msg)))
  }

  def warn(msg: String, t: Throwable): Unit = {
    workflowFileLogger.warn(format(msg), t)
    akkaLogger.foreach(_.warning(format(msg, t)))
    otherLoggers.foreach(_.warn(format(msg), t))
  }

  def error(msg: String): Unit = {
    workflowFileLogger.error(format(msg))
    akkaLogger.foreach(_.error(format(msg)))
    otherLoggers.foreach(_.error(format(msg)))
  }

  def error(msg: String, t: Throwable): Unit = {
    workflowFileLogger.error(format(msg), t)
    akkaLogger.foreach(_.error(t, format(msg)))
    otherLoggers.foreach(_.error(format(msg), t))
  }

  def debug(msg: String): Unit = {
    workflowFileLogger.debug(format(msg))
    akkaLogger.foreach(_.debug(format(msg)))
    otherLoggers.foreach(_.debug(format(msg)))
  }

  def debug(msg: String, t: Throwable): Unit = {
    workflowFileLogger.debug(format(msg), t)
    akkaLogger.foreach(_.debug(format(msg, t)))
    otherLoggers.foreach(_.debug(format(msg), t))
  }

  def trace(msg: String): Unit = {
    workflowFileLogger.trace(format(msg))
    otherLoggers.foreach(_.trace(format(msg)))
  }

  def trace(msg: String, t: Throwable): Unit = {
    workflowFileLogger.trace(format(msg), t)
    otherLoggers.foreach(_.trace(format(msg), t))
  }

  def info(msg: String): Unit = {
    workflowFileLogger.info(format(msg))
    akkaLogger.foreach(_.info(format(msg)))
    otherLoggers.foreach(_.info(format(msg)))
  }

  def info(msg: String, t: Throwable): Unit = {
    workflowFileLogger.info(format(msg), t)
    akkaLogger.foreach(_.info(format(msg, t)))
    otherLoggers.foreach(_.info(format(msg), t))
  }
}
