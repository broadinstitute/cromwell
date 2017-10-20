package wdl.exception

import lenthall.exception.ThrowableAggregation

case class ValidationException(exceptionContext: String, throwables: List[Throwable]) extends ThrowableAggregation