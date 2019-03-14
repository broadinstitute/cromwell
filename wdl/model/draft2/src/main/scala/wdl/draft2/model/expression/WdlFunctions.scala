package wdl.draft2.model.expression

import scala.util.{Failure, Try}

trait WdlFunctions[T] {
  type WdlFunction = Seq[Try[T]] => Try[T]

  /* Returns one of the standard library functions (defined above) by name */
  def getFunction(name: String): WdlFunction = {
    val method = getClass.getMethod(name, classOf[Seq[Try[T]]])
    args => method.invoke(this, args).asInstanceOf[Try[T]]
  }

  /**
    * Extract a single `WomValue` from the specified `Seq`, returning `Failure` if the parameters
    * represent something other than a single `WomValue`.
    */
  def extractSingleArgument(functionName: String, params: Seq[Try[T]]): Try[T] = {
    if (params.length != 1) Failure(new UnsupportedOperationException(s"Expected one argument for $functionName, got ${params.length}"))
    else params.head
  }

  /*
   * Below are methods that can be overridden, if necessary, by engine implementations of the standard library
   * to accommodate particularities of the engine's backend.
   */
  /**
    * Path where to write files created by standard functions (write_*).
   */
  def tempFilePath: String = throw new UnsupportedOperationException("write_* functions are not supported by this implementation")
}
