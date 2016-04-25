package cromwell.backend.impl

import lenthall.exception.AggregatedException

import scala.util.{Failure, Success, Try}

package object local {

  // Duplicate of engine.util
  //TODO lenthall-ify
  object TryUtils {
    private def sequenceIterable[T](tries: Iterable[Try[_]], unbox: () => T, prefixErrorMessage: String) = {
      tries collect { case f: Failure[_] => f } match {
        case failures if failures.nonEmpty => Failure(new AggregatedException(prefixErrorMessage, failures.map(_.failed.get)))
        case _ => Success(unbox())
      }
    }

    def sequence[T](tries: Seq[Try[T]], prefixErrorMessage: String = ""): Try[Seq[T]] = {
      def unbox = tries map {
        _.get
      }
      sequenceIterable(tries, unbox _, prefixErrorMessage)
    }

    def sequenceMap[T, U](tries: Map[T, Try[U]], prefixErrorMessage: String = ""): Try[Map[T, U]] = {
      def unbox = tries mapValues {
        _.get
      }
      sequenceIterable(tries.values, unbox _, prefixErrorMessage)
    }
  }
}
