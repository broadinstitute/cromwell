package cromwell.backend.impl.bcs

import cats.data.Validated._
import lenthall.exception.MessageAggregation
import lenthall.validation.ErrorOr._
import cats.syntax.cartesian._
import cats.syntax.validated._

import scala.util.Try
import scala.util.matching.Regex

object BcsUserData {
  val keyPattern = """[^\s]+"""
  val valuePattern = """[^\s]+"""
  val inputMountPattern: Regex = s"""($keyPattern)\\s+($valuePattern)""".r

  def parse(s: String): Try[BcsUserData] = {
    val validation: ErrorOr[BcsUserData] = s match {
      case inputMountPattern(key, value) =>
        (validateKey(key) |@| validateValue(value)) map { (key, value)
        => new BcsUserData(key, value)
        }
      case _ => s"error user data entry".invalidNel
    }

    Try(validation match {
      case Valid(userData) => userData
      case Invalid(nels) =>
        throw new UnsupportedOperationException with MessageAggregation {
          val exceptionContext = ""
          val errorMessages: List[String] = nels.toList
        }
    })
  }

  private def validateKey(key: String): ErrorOr[String] = {
    key.validNel
  }

  private def validateValue(value: String): ErrorOr[String] = {
    value.validNel
  }

}

case class BcsUserData(key: String, value: String)
