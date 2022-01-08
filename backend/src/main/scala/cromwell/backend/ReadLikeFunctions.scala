package cromwell.backend

import cromwell.core.io.AsyncIoFunctions
import cromwell.core.path.PathFactory
import wom.expression.IoFunctionSet

import scala.concurrent.Future
import scala.util.Try

trait ReadLikeFunctions extends PathFactory with IoFunctionSet with AsyncIoFunctions {

  override def readFile(path: String, maxBytes: Option[Int], failOnOverflow: Boolean): Future[String] =
    Future.fromTry(Try(buildPath(path))) flatMap { p => asyncIo.contentAsStringAsync(p, maxBytes, failOnOverflow) }

  override def size(path: String): Future[Long] = Future.fromTry(Try(buildPath(path))) flatMap asyncIo.sizeAsync
}
