package wom.expression
import scala.concurrent.Future

class IoFunctionSetAdapter(delegate: IoFunctionSet) extends IoFunctionSet {
  override def pathFunctions: PathFunctionSet = delegate.pathFunctions
  override def resolvedPath(path: String): Future[String] = delegate.resolvedPath(path)
  override def readFile(path: String, maxBytes: Option[Int], failOnOverflow: Boolean): Future[String] = delegate.readFile(path, maxBytes, failOnOverflow)
  override def writeFile(path: String, content: String) = delegate.writeFile(path, content)
  override def createTemporaryDirectory(name: Option[String]) = delegate.createTemporaryDirectory(name)
  override def copyFile(source: String, destination: String) = delegate.copyFile(source, destination)
  override def glob(pattern: String) = delegate.glob(pattern)
  override def listAllFilesUnderDirectory(dirPath: String) = delegate.listAllFilesUnderDirectory(dirPath)
  override def listDirectory(path: String)(visited: Vector[String]) = delegate.listDirectory(path)(visited)
  override def isDirectory(path: String) = delegate.isDirectory(path)
  override def size(path: String) = delegate.size(path)
  override implicit def ec  = delegate.ec
}
