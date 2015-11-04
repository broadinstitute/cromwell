package cromwell

import scala.util.Try

/**
 * DSL for processing Main.run() arguments.
 */
object MainRunDsl {

  /*
   * "Should" classes help the compiler anchor the "to" methods, after java jettisons away the generics.  However,
   * erasure + overload still necessitate that the method arg append "_" to become a partially applied function.
   */

  abstract class ShouldOnPath(action: String) {
    override def toString = action
  }

  abstract class ShouldOnJson(action: String) {
    override def toString = action
  }

  case object ShouldParse extends ShouldOnJson("parse")

  case object ShouldProcess extends ShouldOnPath("process")

  case object ShouldAccess extends ShouldOnPath("access")

  /**
   * Try to run an operation, or return a Failure, printing what we were trying to do.
   */
  object Trying {
    /**
     * Try to run a function and return the output.
     * @param functionOn Function to run.
     * @param path Path to run on.
     * @param should A description of what action we're performing.
     * @param inputDescription The description of the input path.
     * @tparam PathType The type of path. Ex: Path or Option[Path].
     * @tparam OutputType The return type of the attempt.
     * @return The Try result of the attempt.
     */
    def to[PathType, OutputType](functionOn: (String, PathType) => OutputType,
                                 path: PathType,
                                 should: ShouldOnPath,
                                 inputDescription: String): Try[OutputType] = {
      attempt(Try(functionOn(inputDescription, path)), s"$should ${inputDescription.toLowerCase}")
    }

    /**
     * Try to run a processing function on some json String.
     * @param tryFunctionOn The function that will attempt to consume the json.
     * @param json The json.
     * @param should A description of what action we're performing. Currently only ShouldParse.
     * @param inputDescription The description of the input json.
     * @return The Try result of the attempt.
     */
    def to(tryFunctionOn: (String) => Try[String],
           json: String,
           should: ShouldOnJson,
           inputDescription: String): Try[String] = {
      attempt(tryFunctionOn(json), s"$should ${inputDescription.toLowerCase} json")
    }

    /**
     * Try to run a function. If it fails, report back to the user.
     * @param tryFunction Function to run.
     * @param stepDescription The description of what the function is attempting.
     * @tparam T The type returned upon success.
     * @return The Try result of the attempt.
     */
    private def attempt[T](tryFunction: => Try[T], stepDescription: String): Try[T] = {
      val attempt = tryFunction
      if (attempt.isFailure) Console.err.println(s"ERROR: Unable to $stepDescription")
      attempt
    }
  }

}
