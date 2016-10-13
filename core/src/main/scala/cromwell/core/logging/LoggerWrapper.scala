package cromwell.core.logging

import akka.event.LoggingAdapter
import org.apache.commons.lang3.exception.ExceptionUtils
import org.slf4j.Logger
import org.slf4j.helpers.MarkerIgnoringBase

/**
  * Wraps an akka logger and a Set of slf4j loggers together.
  */
abstract class LoggerWrapper extends MarkerIgnoringBase {

  def akkaLogger: Option[LoggingAdapter]
  def slf4jLoggers: Set[Logger]

  def tag: String

  private def format(msg: String): String = s"$tag: $msg"

  private def format(msg: String, throwable: Throwable): String = {
    format(msg) + ":\n" + ExceptionUtils.getStackTrace(throwable)
  }

  /**
    * Akka logger only supports 1 to 4 arguments, unlike sl4j logger which support arbitrary number via varargs.
    * This method matches the varargs to the correct akka logger method, if it exists.
    */
  private def varargsAkkaLog(pattern: String, arguments: Any*): Unit = {
    lazy val formatted: String = format(pattern)

    arguments.length match {
      case 1 => akkaLogger.foreach(_.error(formatted, arguments.head))
      case 2 => akkaLogger.foreach(_.error(formatted, arguments.head, arguments(1)))
      case 3 => akkaLogger.foreach(_.error(formatted, arguments.head, arguments(1), arguments(2)))
      case 4 => akkaLogger.foreach(_.error(formatted, arguments.head, arguments(1), arguments(2), arguments(3)))
      case _ =>
    }
  }

  override def warn(pattern: String, arguments: AnyRef*): Unit = {
    lazy val formatted: String = format(pattern)

    varargsAkkaLog(pattern, arguments)
    // The :_* indicates the scala compiler to use the (String, Object...) method and not the (String, Object)
    slf4jLoggers.foreach(_.warn(formatted, arguments:_*))
  }

  override def warn(pattern: String, arg1: Any, arg2: Any): Unit = {
    lazy val formatted: String = format(pattern)

    akkaLogger.foreach(_.warning(formatted, arg1, arg2))
    // The :Any for arg2 indicates the scala compiler to use the (String, Object, Object) method and not the (String, Object...)
    // See https://github.com/typesafehub/scalalogging/issues/16
    slf4jLoggers.foreach(_.warn(formatted, arg1, arg2: Any))
  }

  override def warn(pattern: String, argument: Any): Unit = {
    lazy val formatted: String = format(pattern)

    akkaLogger.foreach(_.warning(formatted, argument))
    slf4jLoggers.foreach(_.warn(formatted, argument))
  }

  override def warn(msg: String): Unit = {
    lazy val formatted: String = format(msg)

    akkaLogger.foreach(_.warning(formatted))
    slf4jLoggers.foreach(_.warn(formatted))
  }

  override def warn(msg: String, t: Throwable): Unit = {
    lazy val formatted: String = format(msg)

    akkaLogger.foreach(_.warning(format(msg, t)))
    slf4jLoggers.foreach(_.warn(formatted, t))
  }

  override def error(msg: String): Unit = {
    lazy val formatted: String = format(msg)

    akkaLogger.foreach(_.error(formatted))
    slf4jLoggers.foreach(_.error(formatted))
  }

  override def error(msg: String, t: Throwable): Unit = {
    lazy val formatted: String = format(msg)

    akkaLogger.foreach(_.error(t, formatted))
    slf4jLoggers.foreach(_.error(formatted, t))
  }

  override def error(pattern: String, arguments: AnyRef*): Unit = {
    lazy val formatted: String = format(pattern)

    varargsAkkaLog(pattern, arguments)
    slf4jLoggers.foreach(_.error(formatted, arguments:_*))
  }

  override def error(pattern: String, arg: Any): Unit = {
    lazy val formatted: String = format(pattern)

    akkaLogger.foreach(_.error(formatted, arg))
    slf4jLoggers.foreach(_.error(formatted, arg))
  }

  override def error(pattern: String, arg1: Any, arg2: Any): Unit = {
    lazy val formatted: String = format(pattern)

    akkaLogger.foreach(_.error(formatted, arg1, arg2))
    slf4jLoggers.foreach(_.error(formatted, arg1, arg2: Any))
  }

  def error(t: Throwable, pattern: String, arguments: Any*): Unit = {
    lazy val formatted: String = format(pattern)

    akkaLogger.foreach(_.error(t, formatted, arguments))
    slf4jLoggers.foreach(_.error(format(pattern, t), arguments))
  }

  override def debug(msg: String): Unit = {
    lazy val formatted: String = format(msg)

    akkaLogger.foreach(_.debug(formatted))
    slf4jLoggers.foreach(_.debug(formatted))
  }

  override def debug(msg: String, t: Throwable): Unit = {
    lazy val formatted: String = format(msg)

    akkaLogger.foreach(_.debug(format(msg, t)))
    slf4jLoggers.foreach(_.debug(formatted, t))
  }

  override def debug(pattern: String, arguments: AnyRef*): Unit = {
    lazy val formatted: String = format(pattern)

    varargsAkkaLog(pattern, arguments)
    slf4jLoggers.foreach(_.debug(formatted, arguments:_*))
  }

  override def debug(pattern: String, argument: Any): Unit = {
    lazy val formatted: String = format(pattern)

    akkaLogger.foreach(_.debug(formatted, argument))
    slf4jLoggers.foreach(_.debug(formatted, argument))
  }

  override def debug(pattern: String, arg1: Any, arg2: Any): Unit = {
    lazy val formatted: String = format(pattern)

    akkaLogger.foreach(_.debug(formatted, arg1, arg2))
    slf4jLoggers.foreach(_.debug(formatted, arg1, arg2: Any))
  }

  override def trace(msg: String): Unit = {
    slf4jLoggers.foreach(_.trace(format(msg)))
  }

  override def trace(msg: String, t: Throwable): Unit = {
    slf4jLoggers.foreach(_.trace(format(msg), t))
  }

  override def trace(pattern: String, arguments: AnyRef*): Unit = {
    slf4jLoggers.foreach(_.trace(format(pattern), arguments:_*))
  }

  override def trace(pattern: String, arg: Any): Unit = {
    slf4jLoggers.foreach(_.trace(format(pattern), arg))
  }

  override def trace(pattern: String, arg1: Any, arg2: Any): Unit = {
    slf4jLoggers.foreach(_.trace(format(pattern), arg1, arg2: Any))
  }

  override def info(msg: String): Unit = {
    lazy val formatted: String = format(msg)

    akkaLogger.foreach(_.info(formatted))
    slf4jLoggers.foreach(_.info(formatted))
  }

  override def info(msg: String, t: Throwable): Unit = {
    lazy val formatted: String = format(msg)

    akkaLogger.foreach(_.info(format(msg, t)))
    slf4jLoggers.foreach(_.info(formatted, t))
  }

  override def info(pattern: String, arguments: AnyRef*): Unit = {
    lazy val formatted: String = format(pattern)

    varargsAkkaLog(pattern, arguments)
    slf4jLoggers.foreach(_.info(formatted, arguments:_*))
  }

  override def info(pattern: String, arg: Any): Unit = {
    lazy val formatted: String = format(pattern)

    akkaLogger.foreach(_.info(formatted, arg))
    slf4jLoggers.foreach(_.info(formatted, arg))
  }

  override def info(pattern: String, arg1: Any, arg2: Any): Unit = {
    lazy val formatted: String = format(pattern)

    akkaLogger.foreach(_.info(formatted, arg1, arg2))
    slf4jLoggers.foreach(_.info(formatted, arg1, arg2: Any))
  }
  
  override def isErrorEnabled: Boolean = throw new UnsupportedOperationException("This logger wraps an arbitrary set of loggers that can each have a different level enabled.")

  override def isInfoEnabled: Boolean = throw new UnsupportedOperationException("This logger wraps an arbitrary set of loggers that can each have a different level enabled.")

  override def isDebugEnabled: Boolean = throw new UnsupportedOperationException("This logger wraps an arbitrary set of loggers that can each have a different level enabled.")

  override def isTraceEnabled: Boolean = throw new UnsupportedOperationException("This logger wraps an arbitrary set of loggers that can each have a different level enabled.")

  override def isWarnEnabled: Boolean = throw new UnsupportedOperationException("This logger wraps an arbitrary set of loggers that can each have a different level enabled.")
  
}
