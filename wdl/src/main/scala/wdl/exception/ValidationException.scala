package wdl.exception

import common.exception.ThrowableAggregation

case class ValidationException(exceptionContext: String, throwables: List[Throwable]) extends ThrowableAggregation