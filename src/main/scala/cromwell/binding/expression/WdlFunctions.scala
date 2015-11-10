package cromwell.binding.expression

import cromwell.binding.values.WdlValue

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

  /**
    * Given a WDL value that represents a file, return the contents
    * of that file.  Not all WdlValues can be interpreted as files
    * In which case this method should throw an exception
    *
    * @param value - WDL Value that represents a file
    * @return - Contents of the file
    * @throws UnsupportedOperationException if the WDL value can
    *         not be interpreted as a file
    * @throws NotImplementedError if the backend did not implement
    *         this method
    */
  def fileContentsToString(value: WdlValue): String = throw new NotImplementedError("fileContentsToString() is unimplemented")

  /* Returns one of the standard library functions (defined above) by name */
  def getFunction(name: String): WdlFunction = {
    val method = getClass.getMethod(name, classOf[Seq[Try[T]]])
    args => method.invoke(this, args).asInstanceOf[Try[T]]
  }
}

