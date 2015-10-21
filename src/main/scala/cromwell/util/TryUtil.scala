package cromwell.util

import java.io.{PrintWriter, StringWriter}

import com.typesafe.scalalogging.LazyLogging

import scala.concurrent.duration.Duration
import scala.util.{Success, Failure, Try}

case class AggregatedException[A](exceptions: Seq[Failure[A]]) extends Exception {
  override def getMessage: String = {
    exceptions.map(_.exception.getMessage).mkString("\n")
  }
}

object TryUtil extends LazyLogging {
  private def stringifyFailure[T](failure: Try[T]): String = {
    val stringWriter = new StringWriter()
    val writer = new PrintWriter(stringWriter)
    failure.recover { case e => e.printStackTrace(writer)}
    writer.flush()
    writer.close()
    stringWriter.toString
  }

  def stringifyFailures[T](possibleFailures: Traversable[Try[T]]): Traversable[String] =
    possibleFailures.collect { case failure: Failure[T] => stringifyFailure(failure) }

  private def defaultSuccessFunction(a: Any): Boolean = true

  /**
   * Runs a block of code (`fn`) `retries` number of times until it succeeds.
   * It will wait `pollingInterval` amount of time between retry attempts and
   * The `pollingBackOffFactor` is for exponentially backing off the `pollingInterval`
   * on subsequent retries.  The `pollingInterval` shall not exceed `maxPollingInterval`
   *
   * Returns a Try[T] where T is the return value of `fn`, the function to be retried.
   * If the return value is Success[T] then at least one retry succeeded.
   *
   * The isSuccess function is optional but if provided, then isSuccess(fn) must be true
   * or it will trigger another retry.  if isSuccess is omitted, the only way the fn can
   * fail is if it throws an exception.
   *
   * Use `retries` value of None indicates to retry indefinitely.
   */
  @annotation.tailrec
  def retryBlock[T](fn: Option[T] => T,
                    isSuccess: T => Boolean = defaultSuccessFunction _,
                    retries: Option[Int],
                    pollingInterval: Duration,
                    pollingBackOffFactor: Double,
                    maxPollingInterval: Duration,
                    failMessage: Option[String] = None,
                    priorValue: Option[T] = None): Try[T] = {
    Try { fn(priorValue) } match {
      case Success(x) if isSuccess(x) => Success(x)
      case value if (retries.isDefined && retries.get > 1) || retries.isEmpty =>

        val retryCountMessage = if (retries.getOrElse(0) > 0) s" (${retries.getOrElse(0) - 1} more retries) " else ""
        val retryMessage = s"Retrying in $pollingInterval$retryCountMessage..."
        failMessage foreach { m => logger.warn(s"$m.  $retryMessage") }

        Thread.sleep(pollingInterval.toMillis)

        retryBlock(
          fn,
          isSuccess,
          retries.map(_ - 1),
          Duration(Math.min((pollingInterval.toMillis * pollingBackOffFactor).toLong, maxPollingInterval.toMillis), "milliseconds"),
          pollingBackOffFactor,
          maxPollingInterval,
          failMessage,
          value.toOption
        )
      case f => f
    }
  }

  def flatten[A](s: Seq[Try[A]]): Try[Seq[A]] = {
    s.collect({case f: Failure[_] => f}) match {
      case failures if failures.nonEmpty => Failure(new AggregatedException(failures))
      case _ => Success(s.map(_.get))
    }
  }
}
