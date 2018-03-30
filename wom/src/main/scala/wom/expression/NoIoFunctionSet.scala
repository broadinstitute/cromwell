package wom.expression

import java.util.concurrent.Executors

import wom.values.WomSingleFile

import scala.concurrent.{ExecutionContext, Future}

class EmptyIoFunctionSet extends IoFunctionSet {
  override def readFile(path: String, maxBytes: Option[Int] = None, failOnOverflow: Boolean = false): Future[String] = Future.failed(new NotImplementedError("readFile is not available here"))

  override def writeFile(path: String, content: String): Future[WomSingleFile] = {
    Future.failed(new NotImplementedError("writeFile is not available here"))
  }

  override def copyFile(pathFrom: String, targetName: String): Future[WomSingleFile] = {
    throw new Exception("copyFile is not available here")
  }

  override def glob(pattern: String): Future[Seq[String]] = throw new NotImplementedError("glob is not available here")

  override def listAllFilesUnderDirectory(dirPath: String): Nothing =
    throw new NotImplementedError("listAllFilesUnderDirectory is not available here")

  override def size(path: String): Future[Long] = Future.failed(new NotImplementedError("size is not available here"))

  override implicit def ec: ExecutionContext = ExecutionContext.fromExecutor(Executors.newFixedThreadPool(1))
  override def pathFunctions = NoPathFunctionSet
  override def listDirectory(path: String) = throw new NotImplementedError("listDirectory is not available here")
  override def isDirectory(path: String) = throw new NotImplementedError("isDirectory is not available here")
}

class EmptyPathFunctionSet extends PathFunctionSet {
  override def sibling(of: String, path: String) = throw new NotImplementedError("sibling is not available here")
  override def isAbsolute(path: String) = false
  override def relativeToHostCallRoot(path: String) = path
  override def name(path: String) = throw new NotImplementedError("name is not available here")
  override def stdout = throw new NotImplementedError("stdout is not available here")
  override def stderr = throw new NotImplementedError("stderr is not available here")
}

case object NoIoFunctionSet extends EmptyIoFunctionSet
case object NoPathFunctionSet extends EmptyPathFunctionSet
