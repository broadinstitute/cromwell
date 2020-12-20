package cromwell.core.logging

import akka.event.Logging.LogLevel
import akka.event.{Logging, LoggingAdapter}
import org.apache.commons.lang3.exception.ExceptionUtils
import org.slf4j.Logger
import org.slf4j.helpers.{MarkerIgnoringBase, MessageFormatter}

/**
  * Wraps an akka logger and a Set of slf4j loggers together in slf4j style.
  */
abstract class LoggerWrapper extends MarkerIgnoringBase {

  def akkaLogger: Option[LoggingAdapter]
  def slf4jLoggers: Set[Logger]

  def tag: String

  private def format(msg: String): String = s"$tag: $msg"

  /**
    * Formats a throwable for akka similar to slf4j's SimpleLogger.
    *
    * https://github.com/qos-ch/slf4j/blob/v_1.7.30/slf4j-simple/src/main/java/org/slf4j/impl/SimpleLogger.java#L293-L295
    */
  private def format(msg: String, throwable: Throwable): String = {
    format(msg) + "\n" + ExceptionUtils.getStackTrace(throwable)
  }

  /**
    * Passes a formatted string to akka similar to slf4j's SimpleLogger
    */
  private def varargsAkkaLog(logLevel: LogLevel, pattern: String, arguments: Seq[Any]): Unit = {
    lazy val formatted: String = format(pattern)

    akkaLogger foreach { akkaLoggerActual =>
      if (akkaLoggerActual.isEnabled(logLevel)) {
        arguments.length match {
          case 0 => akkaLoggerActual.log(logLevel, formatted)
          case _ =>
            // Ignore Akka formatting style. Instead format message the way slf4j formats messages.
            val slf4jTuple = MessageFormatter.arrayFormat(formatted, arguments.map(_.asInstanceOf[AnyRef]).toArray)

            (logLevel, Option(slf4jTuple.getThrowable)) match {
              // Special case: When logging exceptions at an error log level, use akka formatting for passing the cause
              case (Logging.ErrorLevel, Some(cause)) => akkaLoggerActual.error(cause, slf4jTuple.getMessage)
              // If a cause was located by slf4j but we're not logging an exception/error, slf4j will still log it.
              // Do the same here.
              case (_, Some(cause)) => akkaLoggerActual.log(logLevel, format(slf4jTuple.getMessage, cause))
              // Otherwise, just do a normal log message.
              case (_, None) => akkaLoggerActual.log(logLevel, slf4jTuple.getMessage)
            }
        }
      }
    }
  }

  override def warn(pattern: String, arguments: AnyRef*): Unit = {
    lazy val formatted: String = format(pattern)

    varargsAkkaLog(Logging.WarningLevel, pattern, arguments)
    // The :_* indicates the scala compiler to use the (String, Object...) method and not the (String, Object)
    slf4jLoggers.foreach { logger =>
      logger.warn(formatted, arguments: _*)
    }
  }

  override def warn(pattern: String, arg1: Any, arg2: Any): Unit = {
    lazy val formatted: String = format(pattern)

    varargsAkkaLog(Logging.WarningLevel, pattern, Seq(arg1, arg2))
    // The :Any for arg2 indicates the scala compiler to use the (String, Object, Object) method and not the (String, Object...)
    // See https://github.com/typesafehub/scalalogging/issues/16
    slf4jLoggers.foreach(_.warn(formatted, arg1, arg2: Any))
  }

  override def warn(pattern: String, argument: Any): Unit = {
    lazy val formatted: String = format(pattern)

    varargsAkkaLog(Logging.WarningLevel, pattern, Seq(argument))
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

    varargsAkkaLog(Logging.ErrorLevel, pattern, arguments)
    slf4jLoggers.foreach(_.error(formatted, arguments:_*))
  }

  override def error(pattern: String, arg: Any): Unit = {
    lazy val formatted: String = format(pattern)

    varargsAkkaLog(Logging.ErrorLevel, pattern, Seq(arg))
    slf4jLoggers.foreach(_.error(formatted, arg))
  }

  override def error(pattern: String, arg1: Any, arg2: Any): Unit = {
    lazy val formatted: String = format(pattern)

    varargsAkkaLog(Logging.ErrorLevel, pattern, Seq(arg1, arg2))
    slf4jLoggers.foreach(_.error(formatted, arg1, arg2: Any))
  }

  def error(t: Throwable, pattern: String, arguments: Any*): Unit = {
    // slf4j extracts the last variable argument as a throwable.
    error(pattern, (arguments :+ t).map(_.asInstanceOf[AnyRef]): _*)
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

    varargsAkkaLog(Logging.DebugLevel, pattern, arguments)
    slf4jLoggers.foreach(_.debug(formatted, arguments:_*))
  }

  override def debug(pattern: String, argument: Any): Unit = {
    lazy val formatted: String = format(pattern)

    varargsAkkaLog(Logging.DebugLevel, pattern, Seq(argument))
    slf4jLoggers.foreach(_.debug(formatted, argument))
  }

  override def debug(pattern: String, arg1: Any, arg2: Any): Unit = {
    lazy val formatted: String = format(pattern)

    varargsAkkaLog(Logging.DebugLevel, pattern, Seq(arg1, arg2))
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

    varargsAkkaLog(Logging.InfoLevel, pattern, arguments)
    slf4jLoggers.foreach(_.info(formatted, arguments:_*))
  }

  override def info(pattern: String, arg: Any): Unit = {
    lazy val formatted: String = format(pattern)

    varargsAkkaLog(Logging.InfoLevel, pattern, Seq(arg))
    slf4jLoggers.foreach(_.info(formatted, arg))
  }

  override def info(pattern: String, arg1: Any, arg2: Any): Unit = {
    lazy val formatted: String = format(pattern)

    varargsAkkaLog(Logging.InfoLevel, pattern, Seq(arg1, arg2))
    slf4jLoggers.foreach(_.info(formatted, arg1, arg2: Any))
  }

  override def isErrorEnabled: Boolean = throw new UnsupportedOperationException("This logger wraps an arbitrary set of loggers that can each have a different level enabled.")

  override def isInfoEnabled: Boolean = throw new UnsupportedOperationException("This logger wraps an arbitrary set of loggers that can each have a different level enabled.")

  override def isDebugEnabled: Boolean = throw new UnsupportedOperationException("This logger wraps an arbitrary set of loggers that can each have a different level enabled.")

  override def isTraceEnabled: Boolean = throw new UnsupportedOperationException("This logger wraps an arbitrary set of loggers that can each have a different level enabled.")

  override def isWarnEnabled: Boolean = throw new UnsupportedOperationException("This logger wraps an arbitrary set of loggers that can each have a different level enabled.")

}
