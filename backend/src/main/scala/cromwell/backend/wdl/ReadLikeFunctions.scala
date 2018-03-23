package cromwell.backend.wdl

import cromwell.core.io.AsyncIoFunctions
import cromwell.core.path.PathFactory
import wdl4s.parser.MemoryUnit
import wom.expression.IoFunctionSet
import wom.format.MemorySize
import wom.types._
import wom.values._

import scala.concurrent.Future
import scala.util.{Success, Try}

trait ReadLikeFunctions extends PathFactory with IoFunctionSet with AsyncIoFunctions {
  
  // TODO WOM: https://github.com/broadinstitute/cromwell/issues/2611
  val fileSizeLimitationConfig = FileSizeLimitationConfig.fileSizeLimitationConfig

  override def readFile(path: String, maxBytes: Option[Int], failOnOverflow: Boolean): Future[String] = asyncIo.contentAsStringAsync(buildPath(path), maxBytes, failOnOverflow)

  override def size(path: String): Future[Long] = asyncIo.sizeAsync(buildPath(path))
}
