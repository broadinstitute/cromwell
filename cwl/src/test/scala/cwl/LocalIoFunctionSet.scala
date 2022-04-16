package cwl

import java.nio.file.{Files, Paths}
import java.util.concurrent.Executors

import wom.expression.EmptyIoFunctionSet

import scala.jdk.CollectionConverters._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

/**
  * Test io functions that can read from local file systems.
  */
object LocalIoFunctionSet extends EmptyIoFunctionSet {
  override implicit def ec: ExecutionContext = ExecutionContext.fromExecutor(Executors.newFixedThreadPool(1))

  private def stripLocalPrefix(path: String): String = {
    path.stripPrefix("file://")
  }

  private def fileSize(path: String): Long = {
    Files.size(Paths.get(stripLocalPrefix(path)))
  }

  private def contentsAsString(path: String): String = {
    Files.readAllLines(Paths.get(stripLocalPrefix(path))).asScala.mkString("\n")
  }

  override def readFile(path: String, maxBytesOption: Option[Int] = None,
                        failOnOverflow: Boolean = false): Future[String] = {
    Future.fromTry(Try {
      maxBytesOption match {
        case Some(maxBytes) if failOnOverflow && fileSize(path) > maxBytes =>
          throw new RuntimeException(s"overflow while reading $path: ${fileSize(path)} > $maxBytes")
        case Some(maxBytes) => contentsAsString(path).take(maxBytes)
        case None => contentsAsString(path)
      }
    })
  }

  override def size(params: String): Future[Long] = {
    Future.fromTry(Try {
        fileSize(params)
    }
    )
  }
}
