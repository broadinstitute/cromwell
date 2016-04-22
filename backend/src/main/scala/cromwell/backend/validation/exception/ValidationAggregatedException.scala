package cromwell.backend.validation.exception

import lenthall.exception.MessageAggregation

case class ValidationAggregatedException(override val exceptionContext: String,
                                         override val errorMessages: Traversable[String]) extends MessageAggregation
