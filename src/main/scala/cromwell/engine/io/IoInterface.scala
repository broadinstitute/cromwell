package cromwell.engine.io

/**
  * Defines the set of IO methods needed by the WdlStandardLibrary to execute IO engine functions (read_*, write_*, glob, ...).
  * Implementation is Backend specific. Authentication might be workflow specific.
  * An instance of this is expected to be created when a WorkflowDescriptor is instantiated.
  * It will then hold a reference to an IOInterface that can be used throughout the workflow execution to evaluate expressions with IO engine functions.
  */
trait IoInterface {
  def readFile(path: String): String
  def writeFile(path: String, content: String): Unit
  def writeTempFile(path: String, prefix: String, suffix: String, content: String): String
  def exists(path: String): Boolean
  def listContents(path: String): Iterable[String]
  def glob(path: String, pattern: String): Seq[String]
  def copy(from: String, to: String): Unit
  def hash(path: String): String
  def isValidPath(path: String): Boolean
}
