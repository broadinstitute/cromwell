package common.validation

import java.net.{URI, URL}

import cats.data.Validated.{Invalid, Valid}
import cats.data.{NonEmptyList, Validated}
import cats.syntax.validated._
import cats.syntax.either._
import common.Checked
import common.exception.AggregatedMessageException
import common.validation.ErrorOr.ErrorOr
import net.ceedubs.ficus.readers.{StringReader, ValueReader}
import org.slf4j.Logger

import scala.util.{Failure, Success, Try}

object Validation {
  def warnNotRecognized(keys: Set[String], reference: Set[String], context: String, logger: Logger): Unit = {
    val unrecognizedKeys = keys.diff(reference)
    if (unrecognizedKeys.nonEmpty) {
      logger.warn(s"Unrecognized configuration key(s) for $context: ${unrecognizedKeys.mkString(", ")}")
    }
  }
  
  implicit val urlReader: ValueReader[URL] = StringReader.stringValueReader.map { URI.create(_).toURL }
  
  def validate[A](block: => A): ErrorOr[A] = Try(block) match {
    case Success(result) => result.validNel
    case Failure(f) => f.getMessage.invalidNel
  }

  implicit class TryValidation[A](val t: Try[A]) extends AnyVal {
    def toErrorOr: ErrorOr[A] = {
      Validated.fromTry(t).leftMap(_.getMessage).toValidatedNel[String, A] 
    }

    def toChecked: Checked[A] = {
      Either.fromTry(t).leftMap(ex => NonEmptyList.of(ex.getMessage))
    }
  }

  implicit class ValidationTry[A](val e: ErrorOr[A]) extends AnyVal {
    def toTry: Try[A] = toTry("Error(s)")

    def toTry(context: String): Try[A] = e match {
      case Valid(options) => Success(options)
      case Invalid(err) => Failure(AggregatedMessageException(context, err.toList))
    }
  }
}
