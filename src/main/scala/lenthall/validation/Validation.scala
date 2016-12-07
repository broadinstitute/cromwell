package lenthall.validation

import java.net.{URI, URL}

import cats.syntax.validated._
import lenthall.validation.ErrorOr.ErrorOr
import net.ceedubs.ficus.readers.{StringReader, ValueReader}
import org.slf4j.Logger

import scala.util.{Failure, Success, Try}

object Validation {
  def warnNotRecognized(keys: Set[String], reference: Set[String], context: String, logger: Logger) = {
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
}
