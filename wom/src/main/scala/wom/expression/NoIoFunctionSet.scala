package wom.expression

import java.util.concurrent.Executors

import wom.values.WomSingleFile

import scala.concurrent.{ExecutionContext, Future}

object EmptyIoFunctionSet {
  lazy val singleThreadEc = ExecutionContext.fromExecutor(Executors.newSingleThreadExecutor())
}

class EmptyIoFunctionSet extends IoFunctionSet {
  override def readFile(path: String, maxBytes: Option[Int] = None, failOnOverflow: Boolean = false): Future[String] = Future.failed(new UnsupportedOperationException("readFile is not available here"))

  override def writeFile(path: String, content: String): Future[WomSingleFile] = {
    Future.failed(new UnsupportedOperationException("writeFile is not available here"))
  }

  override def copyFile(pathFrom: String, targetName: String): Future[WomSingleFile] = {
    throw new Exception("copyFile is not available here")
  }

  override def glob(pattern: String): Future[Seq[String]] = throw new UnsupportedOperationException("glob is not available here")

  override def listAllFilesUnderDirectory(dirPath: String): Nothing =
    throw new UnsupportedOperationException("listAllFilesUnderDirectory is not available here")

  override def size(path: String): Future[Long] = Future.failed(new UnsupportedOperationException("size is not available here"))

  override implicit def ec: ExecutionContext = EmptyIoFunctionSet.singleThreadEc
  override def pathFunctions = NoPathFunctionSet
  override def listDirectory(path: String)(visited: Vector[String] = Vector.empty) = throw new UnsupportedOperationException("listDirectory is not available here")
  override def isDirectory(path: String): Future[Boolean] = throw new UnsupportedOperationException("isDirectory is not available here")
  override def createTemporaryDirectory(name: Option[String]) = throw new UnsupportedOperationException("createTemporaryDirectory is not available here")
}

class EmptyPathFunctionSet extends PathFunctionSet {
  override def sibling(of: String, path: String) = throw new UnsupportedOperationException("sibling is not available here")
  override def isAbsolute(path: String) = false
  override def relativeToHostCallRoot(path: String) = path
  override def name(path: String) = throw new UnsupportedOperationException("name is not available here")
  override def stdout = throw new UnsupportedOperationException("stdout is not available here")
  override def stderr = throw new UnsupportedOperationException("stderr is not available here")
}

case object NoIoFunctionSet extends EmptyIoFunctionSet
case object NoPathFunctionSet extends EmptyPathFunctionSet
