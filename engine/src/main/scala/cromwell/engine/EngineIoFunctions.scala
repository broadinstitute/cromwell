package cromwell.engine

import cromwell.backend.ReadLikeFunctions
import cromwell.core.io.{AsyncIo, WorkflowCorePathFunctions}
import cromwell.core.path.PathBuilder
import wom.values.WomSingleFile
import better.files.File._

import scala.concurrent.{ExecutionContext, Future}

class EngineIoFunctions(val pathBuilders: List[PathBuilder],
                        override val asyncIo: AsyncIo,
                        override val ec: ExecutionContext
) extends ReadLikeFunctions
    with WorkflowCorePathFunctions {
  override def glob(pattern: String): Future[Seq[String]] = throw new UnsupportedOperationException(
    s"glob(path, pattern) not implemented yet"
  )

  // TODO: This is not suited for multi backend / multi filesystem use. Keep local for now to not break local CWL conf tests
  override def writeFile(path: String, content: String): Future[WomSingleFile] = Future.successful {
    val cromwellPath = buildPath(path)
    val string =
      if (cromwellPath.isAbsolute)
        cromwellPath.write(content).pathAsString
      else
        (newTemporaryDirectory() / path).write(content).pathAsString
    WomSingleFile(string)
  }
}
