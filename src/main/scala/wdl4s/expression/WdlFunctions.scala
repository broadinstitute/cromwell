package wdl4s.expression

import scala.language.postfixOps
import scala.util.{Failure, Try}

trait WdlFunctions[T] {
  type WdlFunction = Seq[Try[T]] => Try[T]

  /**
   * Extract a single `WdlValue` from the specified `Seq`, returning `Failure` if the parameters
   * represent something other than a single `WdlValue`.
   */
  protected def extractSingleArgument(params: Seq[Try[T]]): Try[T] = {
    if (params.length != 1) Failure(new UnsupportedOperationException("Expected one argument, got " + params.length))
    else params.head
  }

  /* Returns one of the standard library functions (defined above) by name */
  def getFunction(name: String): WdlFunction = {
    val method = getClass.getMethod(name, classOf[Seq[Try[T]]])
    args => method.invoke(this, args).asInstanceOf[Try[T]]
  }

  /*
   * Below are methods that can be overridden, if necessary, by engine implementations of the standard library
   * to accommodate particularities of the engine's backend.
   */
  /**
    * Path where to write files created by standard functions (write_*).
   */
  def tempFilePath: String = throw new NotImplementedError("write_* functions are not supported by this implementation")

  /**
    * Path where to glob from when the glob standard function evaluates.
   */
  def globPath(glob: String): String = throw new NotImplementedError("glob function is not supported by this implementation")
}

