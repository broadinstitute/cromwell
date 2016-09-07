package lenthall.config

import lenthall.exception.MessageAggregation

import scalaz.NonEmptyList

case class ConfigValidationException(configurationContext: String, errors: NonEmptyList[String]) extends RuntimeException with MessageAggregation {
  def this(context: String, errorMessage: String) = this(context, NonEmptyList(errorMessage))

  override val errorMessages = errors.list.toList
  def exceptionContext = s"Invalid $configurationContext configuration"
}
