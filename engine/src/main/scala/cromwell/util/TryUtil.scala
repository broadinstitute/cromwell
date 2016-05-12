package cromwell.util

import java.io.{PrintWriter, StringWriter}

import cromwell.backend.Backoff
import cromwell.core.{CromwellAggregatedException, CromwellFatalException}
import cromwell.logging.WorkflowLogger

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}


@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
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
  private def defaultIsFatal(t: Throwable): Boolean = false
  private def defaultIsTransient(t: Throwable): Boolean = false


  /**
    * Runs a block of code (`fn`) `retries` number of times until it succeeds.
    * It will use an InitializedBackoff instance to handle backoff timings and wait
    * between each retry.
    *
    * Returns a Try[T] where T is the return value of `fn`, the function to be retried.
    * If the return value is Success[T] then at least one retry succeeded.
    *
    * The isSuccess function is optional but if provided, then isSuccess(fn) must be true
    * or it will trigger another retry.  if isSuccess is omitted, the only way the fn can
    * fail is if it throws an exception.
    *
    * The isFatal function is optional but if provided, then a failure f for which
    * isFatal(f) returns true would terminate the retry loop, even if the retry limit has not been reached.
    *
    * The isTransient function is optional but if provided, then a failure f for which
    * isTransient(f) returns true would not count against the retry limit.
    *
    * Note that if the retry limit is reached, the last failure is wrapped into a CromwellFatalException.
    *
    * Use `retryLimit` value of None indicates to retry indefinitely.
    */
  @annotation.tailrec
  def retryBlock[T](fn: Option[T] => T,
                    backoff: Backoff,
                    isSuccess: T => Boolean = defaultSuccessFunction _,
                    retryLimit: Option[Int],
                    logger: WorkflowLogger,
                    failMessage: Option[String] = None,
                    priorValue: Option[T] = None,
                    isFatal: (Throwable) => Boolean = defaultIsFatal _,
                    isTransient: (Throwable) => Boolean = defaultIsTransient _): Try[T] = {
    def logFailures(attempt: Try[T]): Unit = {
      attempt recover {
        case t: Throwable => logger.warn(t.getMessage, t)
      }
    }

    Try { fn(priorValue) } match {
      case Success(x) if isSuccess(x) => Success(x)
      case Failure(f) if isFatal(f) => Failure(new CromwellFatalException(f))
      case value if (retryLimit.isDefined && retryLimit.get > 1) || retryLimit.isEmpty =>
        logFailures(value)

        val transient = isTransient(value.failed.get)
        val newLimit = if (transient) retryLimit else retryLimit.map(_ - 1)
        val retryCountMessage = {
          val transientMessage = if (transient) " - This error is transient and does not count against the retry limit." else ""
          if (retryLimit.getOrElse(0) > 0 || transient) s" (${newLimit.getOrElse(0)} more retries)$transientMessage" else ""
        }
        val pollingInterval = backoff.backoffMillis
        val retryMessage = s"Retrying in $pollingInterval$retryCountMessage..."
        failMessage foreach { m => logger.warn(s"$m.  $retryMessage") }

        Thread.sleep(pollingInterval)

        retryBlock(
          fn,
          backoff.next,
          isSuccess,
          newLimit,
          logger,
          failMessage,
          value.toOption,
          isFatal = isFatal,
          isTransient = isTransient
        )
      case f =>
        logFailures(f)
        f recoverWith {
          case e => Failure(new CromwellFatalException(e))
        }
    }
  }

  private def sequenceIterable[T](tries: Iterable[Try[_]], unbox: () => T, prefixErrorMessage: String) = {
    tries collect { case f: Failure[_] => f } match {
      case failures if failures.nonEmpty => Failure(new CromwellAggregatedException(failures map { _.exception } toSeq, prefixErrorMessage))
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
