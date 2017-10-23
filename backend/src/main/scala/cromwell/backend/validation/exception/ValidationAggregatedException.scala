package cromwell.backend.validation.exception

import common.exception.MessageAggregation

case class ValidationAggregatedException(override val exceptionContext: String,
                                         override val errorMessages: Traversable[String]) extends MessageAggregation
