package lenthall.exception

import java.io.FileNotFoundException
import java.nio.file.NoSuchFileException
import lenthall.exception.Aggregation._

object Aggregation {
  def formatMessageWithList(message: String, list: Traversable[String]) = {
    if (list.nonEmpty) {
      val messages = s"\n${list.mkString("\n")}"
      s"$message:$messages"
    } else message
  }
}

/**
  * When mixed in an Exception class,
  * aggregates multiple error messages into the getMessage method.
  */
trait MessageAggregation extends Throwable {
  def exceptionContext: String
  def errorMessages: Traversable[String]

  override def getMessage = formatMessageWithList(exceptionContext, errorMessages)
}

/**
  * When mixed in an Exception class,
  * aggregates multiple throwables into the extended Exception.
  */
trait ThrowableAggregation extends MessageAggregation {
  def throwables: Traversable[Throwable]

  throwables foreach addSuppressed

  override def errorMessages = throwables map buildMessage

  private def buildMessage(t: Throwable): String = t match {
      // The message for file not found exception only contains the file name, so add the actual reason
    case _: FileNotFoundException | _: NoSuchFileException => s"File not found ${t.getMessage}"
    case aggregation: ThrowableAggregation =>
      formatMessageWithList(aggregation.exceptionContext, aggregation.throwables.map(buildMessage).map("\t" + _))
    case other =>
      val cause = Option(other.getCause) map { c => s"\n\t${buildMessage(c)}" }  getOrElse ""
      s"${other.getMessage}$cause"
  }
}

/**
  * Generic convenience case class for aggregated exceptions.
  */
case class AggregatedException(exceptionContext: String, throwables: Traversable[Throwable]) extends Exception with ThrowableAggregation
case class AggregatedMessageException(exceptionContext: String, errorMessages: Traversable[String]) extends Exception with MessageAggregation
