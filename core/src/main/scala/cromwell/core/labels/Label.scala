package cromwell.core.labels

import lenthall.validation.ErrorOr.ErrorOr
import cats.data.Validated._
import cats.syntax.cartesian._
import cats.syntax.validated._

import scala.util.matching.Regex

sealed abstract case class Label(key: String, value: String)

object Label {

  val MaxLabelLength = 63
  val LabelKeyRegex = "[a-z]([-a-z0-9]*[a-z0-9])?"
  val LabelValueRegex = "([a-z0-9]*[-a-z0-9]*[a-z0-9])?"

  val LabelExpectationsMessage =
    s"A Label key must match the pattern `$LabelKeyRegex` and a label value must match the pattern `$LabelValueRegex`."

  def validateLabelRegex(s: String, regexAllowed: Regex): ErrorOr[String] = {
    (regexAllowed.pattern.matcher(s).matches, s.length <= MaxLabelLength) match {
      case (true, true) => s.validNel
      case (false, false) => s"Invalid label: `$s` did not match regex $regexAllowed and it is ${s.length} characters. The maximum is $MaxLabelLength.".invalidNel
      case (false, _) => s"Invalid label: `$s` did not match the regex $regexAllowed.".invalidNel
      case (_, false) => s"Invalid label: `$s` is ${s.length} characters. The maximum is $MaxLabelLength.".invalidNel
    }
  }

  def validateLabelKey(key: String): ErrorOr[String] = validateLabelRegex(key, LabelKeyRegex.r)

  def validateLabelValue(key: String): ErrorOr[String] = validateLabelRegex(key, LabelValueRegex.r)

  def validateLabel(key: String, value: String): ErrorOr[Label] = {
    val validatedKey = validateLabelKey(key)
    val validatedValue = validateLabelValue(value)

    (validatedKey |@| validatedValue) map { case (k, v) => new Label(k, v) {} }
  }

  def apply(key: String, value: String) = {
    new Label(key, value) {}
  }
}
