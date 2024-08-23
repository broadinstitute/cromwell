package wom.expression

class IoFunctionSetAdapter(delegate: IoFunctionSet) extends IoFunctionSet {
  override def pathFunctions = delegate.pathFunctions
  override def readFile(path: String, maxBytes: Option[Int], failOnOverflow: Boolean) =
    delegate.readFile(path, maxBytes, failOnOverflow)
  override def writeFile(path: String, content: String) = delegate.writeFile(path, content)
  override def glob(pattern: String) = delegate.glob(pattern)
  override def size(path: String) = delegate.size(path)
  implicit override def ec = delegate.ec
}
