package cromwell.core.exception

import lenthall.exception.ThrowableAggregation

case class CromwellAggregatedException(throwables: Seq[Throwable], exceptionContext: String = "") extends ThrowableAggregation
