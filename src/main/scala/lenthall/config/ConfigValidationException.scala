package lenthall.config

import lenthall.exception.MessageAggregation

import cats.data.NonEmptyList

case class ConfigValidationException(configurationContext: String, errors: NonEmptyList[String]) extends RuntimeException with MessageAggregation {
  def this(context: String, errorMessage: String) = this(context, NonEmptyList.of(errorMessage))

  override val errorMessages = errors.toList
  def exceptionContext = s"Invalid $configurationContext configuration"
}
