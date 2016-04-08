package lenthall.exception

/**
  * When mixed in an Exception class,
  * aggregates multiple error messages into the getMessage method.
  */
trait MessageAggregation extends Throwable {
  def exceptionContext: String
  def errorMessages: Traversable[String]

  override def getMessage = {
    s"""$exceptionContext\n${errorMessages.mkString("\n")}"""
  }
}

/**
  * When mixed in an Exception class,
  * aggregates multiple throwables into the extended Exception.
  */
trait ThrowableAggregation extends MessageAggregation {
  def throwables: Traversable[Throwable]

  throwables foreach addSuppressed

  override val errorMessages = throwables map { _.getMessage }
}

/**
  * Generic convenience case class for aggregated exceptions.
  */
case class AggregatedException(exceptionContext: String, throwables: Traversable[Throwable]) extends Exception with ThrowableAggregation
