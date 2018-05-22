package cromwell.backend.wdl

import cromwell.core.io.AsyncIoFunctions
import cromwell.core.path.PathFactory
import wom.expression.IoFunctionSet
import scala.concurrent.Future

trait ReadLikeFunctions extends PathFactory with IoFunctionSet with AsyncIoFunctions {
  
  // TODO WOM: https://github.com/broadinstitute/cromwell/issues/2611
  val fileSizeLimitationConfig = FileSizeLimitationConfig.fileSizeLimitationConfig

  override def readFile(path: String, maxBytes: Option[Int], failOnOverflow: Boolean): Future[String] = asyncIo.contentAsStringAsync(buildPath(path), maxBytes, failOnOverflow)

  override def size(path: String): Future[Long] = asyncIo.sizeAsync(buildPath(path))
}
