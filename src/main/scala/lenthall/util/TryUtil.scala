package lenthall.util

import java.io.{PrintWriter, StringWriter}

import lenthall.exception.AggregatedException

import scala.util.{Failure, Success, Try}

object TryUtil {
  private def stringifyFailure[T](failure: Try[T]): String = {
    val stringWriter = new StringWriter()
    val writer = new PrintWriter(stringWriter)
    failure recover { case e => e.printStackTrace(writer) }
    writer.flush()
    writer.close()
    stringWriter.toString
  }

  def stringifyFailures[T](possibleFailures: Traversable[Try[T]]): Traversable[String] =
    possibleFailures.collect { case failure: Failure[T] => stringifyFailure(failure) }

  private def sequenceIterable[T](tries: Iterable[Try[_]], unbox: () => T, prefixErrorMessage: String) = {
    tries collect { case f: Failure[_] => f } match {
      case failures if failures.nonEmpty =>
        val exceptions = failures.toSeq.map(_.exception)
        Failure(AggregatedException(prefixErrorMessage, exceptions.toList))
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

