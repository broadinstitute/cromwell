package cromwell.core.labels

import common.validation.ErrorOr.ErrorOr
import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._

sealed abstract case class Label(key: String, value: String)

object Label {

  val MaxLabelLength = 255
  val LabelExpectationsMessage = s"A Label key must be non-empty."

  def validateLabelKey(s: String): ErrorOr[String] =
    (s.length >= 1, s.length <= MaxLabelLength) match {
      case (true, true) => s.validNel
      case (false, _) => s"Invalid label: `$s` can't be empty".invalidNel
      case (_, false) => s"Invalid label: `$s` is ${s.length} characters. The maximum is $MaxLabelLength.".invalidNel
    }

  def validateLabelValue(s: String): ErrorOr[String] =
    if (s.length <= MaxLabelLength) {
      s.validNel
    } else {
      s"Invalid label: `$s` is ${s.length} characters. The maximum is $MaxLabelLength.".invalidNel
    }

  def validateLabel(key: String, value: String): ErrorOr[Label] = {
    val validatedKey = validateLabelKey(key)
    val validatedValue = validateLabelValue(value)

    (validatedKey, validatedValue) mapN Label.apply
  }

  def apply(key: String, value: String) =
    new Label(key, value) {}
}
