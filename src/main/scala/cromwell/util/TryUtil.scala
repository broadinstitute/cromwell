package cromwell.util

import java.io.{PrintWriter, StringWriter}

import cromwell.logging.WorkflowLogger

import scala.concurrent.duration.{Duration, _}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

case class AggregatedException(exceptions: Seq[Throwable], prefixError: String = "") extends Exception {
  override def getMessage: String = {
    prefixError + exceptions.map(_.getMessage).mkString("\n")
  }
}

object TryUtil {
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
   * Use `retryLimit` value of None indicates to retry indefinitely.
   */
  @annotation.tailrec
  def retryBlock[T](fn: Option[T] => T,
                    isSuccess: T => Boolean = defaultSuccessFunction _,
                    retryLimit: Option[Int],
                    pollingInterval: Duration,
                    pollingBackOffFactor: Double,
                    maxPollingInterval: Duration,
                    logger: WorkflowLogger,
                    failMessage: Option[String] = None,
                    priorValue: Option[T] = None,
                    fatalExceptions: Seq[Class[_ <: Throwable]] = Seq.empty): Try[T] = {

    def logFailures(attempt: Try[T]): Unit = {
      attempt recover {
        case t: Throwable => logger.warn(t.getMessage, t)
      }
    }

    Try { fn(priorValue) } match {
      case Success(x) if isSuccess(x) => Success(x)
      case Failure(f) if fatalExceptions.contains(f.getClass) => Failure(f)
      case value if (retryLimit.isDefined && retryLimit.get > 1) || retryLimit.isEmpty =>
        logFailures(value)
        val retryCountMessage = if (retryLimit.getOrElse(0) > 0) s" (${retryLimit.getOrElse(0) - 1} more retries) " else ""
        val retryMessage = s"Retrying in $pollingInterval$retryCountMessage..."
        failMessage foreach { m => logger.warn(s"$m.  $retryMessage") }

        Thread.sleep(pollingInterval.toMillis)

        retryBlock(
          fn,
          isSuccess,
          retryLimit.map(_ - 1),
          Duration(Math.min((pollingInterval.toMillis * pollingBackOffFactor).toLong, maxPollingInterval.toMillis), "milliseconds"),
          pollingBackOffFactor,
          maxPollingInterval,
          logger,
          failMessage,
          value.toOption
        )
      case f =>
        logFailures(f)
        f
    }
  }

  private def sequenceIterable[T](tries: Iterable[Try[_]], unbox: () => T, prefixErrorMessage: String) = {
    tries collect { case f: Failure[_] => f } match {
      case failures if failures.nonEmpty => Failure(new AggregatedException(failures map { _.exception } toSeq, prefixErrorMessage))
      case _ => Success(unbox())
    }
  }

  def sequence[T](tries: Seq[Try[T]], prefixErrorMessage: String = ""): Try[Seq[T]] = {
    def unbox = tries map { _.get }
    sequenceIterable(tries, unbox _, prefixErrorMessage)
  }

  def sequenceMap[T, U](tries: Map[T, Try[U]], prefixErrorMessage: String = ""): Try[Map[T, U]] = {
    def unbox = tries mapValues { _.get }
    sequenceIterable(tries.values, unbox _, prefixErrorMessage)
  }

}
