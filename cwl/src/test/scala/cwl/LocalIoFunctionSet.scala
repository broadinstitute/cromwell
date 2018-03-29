package cwl

import java.nio.file.{Files, Paths}

import wom.expression.IoFunctionSet
import wom.values.{WomFloat, WomSingleFile, WomValue}

import scala.collection.JavaConverters._
import scala.concurrent.Future
import scala.util.{Success, Try}

/**
  * Test io functions that can read from local file systems.
  */
object LocalIoFunctionSet extends IoFunctionSet {
  private def unsupported(description: String): Nothing = {
    throw new UnsupportedOperationException(description)
  }

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

  override def writeFile(path: String, content: String): Future[WomSingleFile] = unsupported(s"write file $path")

  override def copyFile(pathFrom: String, targetName: String): Future[WomSingleFile] = {
    unsupported(s"copy file $pathFrom to $targetName")
  }

  override def stdout(params: Seq[Try[WomValue]]) = unsupported(s"stdout for $params")

  override def stderr(params: Seq[Try[WomValue]]) = unsupported(s"stderr for $params")

  override def glob(pattern: String): Future[Seq[String]] = unsupported(s"glob for $pattern")

  override def listAllFilesUnderDirectory(dirPath: String): Future[Seq[String]] = {
    unsupported(s"listAllFilesUnderDirectory $dirPath")
  }

  override def size(params: String): Future[Long] = {
    Future.fromTry(Try {
        fileSize(params)
    }
    )
  }
}
