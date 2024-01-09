package common.util

import java.io.{PrintWriter, StringWriter}

import common.exception.AggregatedException
import common.collections.EnhancedCollections._

import scala.util.{Failure, Success, Try}

object TryUtil {
  private def stringifyFailure[T](failure: Failure[T]): String = {
    val stringWriter = new StringWriter()
    val writer = new PrintWriter(stringWriter)
    failure match {
      case Failure(e) => e.printStackTrace(writer)
    }
    writer.flush()
    writer.close()
    stringWriter.toString
  }

  def stringifyFailures[T](possibleFailures: Iterable[Try[T]]): Iterable[String] =
    possibleFailures.collect { case failure: Failure[T] => stringifyFailure(failure) }

  private def sequenceIterable[T](tries: Iterable[Try[_]], unbox: () => T, prefixErrorMessage: String): Try[T] =
    tries collect { case f: Failure[_] => f } match {
      case failures if failures.nonEmpty =>
        val exceptions = failures.toSeq.map(_.exception)
        Failure(AggregatedException(prefixErrorMessage, exceptions.toList))
      case _ => Success(unbox())
    }

  def sequence[T](tries: Seq[Try[T]], prefixErrorMessage: String = ""): Try[Seq[T]] = {
    def unbox = tries map { _.get }
    sequenceIterable(tries, () => unbox, prefixErrorMessage)
  }

  def sequenceOption[T](tried: Option[Try[T]], prefixErrorMessage: String = ""): Try[Option[T]] = {
    def unbox = tried.map(_.get)

    sequenceIterable(tried.toSeq, () => unbox, prefixErrorMessage)
  }

  def sequenceMap[T, U](tries: Map[T, Try[U]], prefixErrorMessage: String = ""): Try[Map[T, U]] = {
    def unbox = tries safeMapValues { _.get }
    sequenceIterable(tries.values, () => unbox, prefixErrorMessage)
  }

  // NOTE: Map is invariant on the key type, so we accept _ <: Try[T]
  def sequenceKeyValues[T, U](tries: Map[_ <: Try[T], Try[U]], prefixErrorMessage: String = ""): Try[Map[T, U]] = {
    def unbox: Map[T, U] = tries map { case (tryKey, tryValue) => tryKey.get -> tryValue.get }

    sequenceIterable(tries.toSeq.flatMap(Function.tupled(Seq(_, _))), () => unbox, prefixErrorMessage)
  }

  def sequenceTuple[T, U](tries: (Try[T], Try[U]), prefixErrorMessage: String = ""): Try[(T, U)] = {
    def unbox: (T, U) = tries match {
      case (try1, try2) => (try1.get, try2.get)
    }

    sequenceIterable(tries match { case (try1, try2) => Seq(try1, try2) }, () => unbox, prefixErrorMessage)
  }
}
