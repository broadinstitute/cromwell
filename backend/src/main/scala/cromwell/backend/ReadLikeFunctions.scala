package cromwell.backend

import cromwell.core.io.AsyncIoFunctions
import cromwell.core.path.PathFactory
import wom.expression.IoFunctionSet

import scala.concurrent.Future
import scala.util.{Failure, Success, Try}

trait ReadLikeFunctions extends PathFactory with IoFunctionSet with AsyncIoFunctions {
  
  override def readFile(path: String, maxBytes: Option[Int], failOnOverflow: Boolean): Future[String] = {
    println("about to read file")

    Future.fromTry(Try(buildPath(path))) flatMap { p =>
      val result = asyncIo.contentAsStringAsync(p, maxBytes, failOnOverflow)
      result.onComplete {
        case Success(value) => println(s"Successfully read: $value")
        case Failure(value) => println(s"Failed with: $value")
      }
      result
    }
  }

  override def size(path: String): Future[Long] = Future.fromTry(Try(buildPath(path))) flatMap asyncIo.sizeAsync
}
