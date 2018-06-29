package cromwell.languages.util

import java.nio.file.Paths

import better.files.File
import cats.data.NonEmptyList
import cats.effect.IO
import cats.syntax.either._
import cats.syntax.validated._
import com.softwaremill.sttp._
import com.softwaremill.sttp.asynchttpclient.cats.AsyncHttpClientCatsBackend
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.core.path.Path
import wom.core.WorkflowSource

import scala.concurrent.duration._
import scala.concurrent.Await
import scala.util.Try

object ImportResolver {
  type ImportResolver = CheckedAtoB[String, WorkflowSource]

  def zippedImportsResolver(zippedImports: Array[Byte]): ImportResolver = {
    directoryResolver(LanguageFactoryUtil.validateImportsDirectory(zippedImports).toEither)
  }

  def directoryResolver(directory: Path, allowEscapingDirectory: Boolean = false): ImportResolver = CheckedAtoB.fromErrorOr { path =>
    Try(Paths.get(directory.resolve(path).toFile.getCanonicalPath)).toErrorOr flatMap { absolutePathToFile =>
      val absolutePathToImports = Paths.get(directory.toJava.getCanonicalPath)
      if (allowEscapingDirectory || absolutePathToFile.startsWith(absolutePathToImports)) {
        val file = File(absolutePathToFile)
        if (file.exists) {
          File(absolutePathToFile).contentAsString.validNel
        } else {
          s"Import file not found: $path".invalidNel
        }
      } else {
        s"$path is not a valid import".invalidNel
      }
    }
  }

  def directoryResolver(directoryValidation: Checked[Path]): ImportResolver = CheckedAtoB.fromCheck { path =>
    directoryValidation flatMap { directoryResolver(_).run(path) }
  }

  lazy val localFileResolver: ImportResolver = CheckedAtoB.fromErrorOr { path =>
    Try(Paths.get(Paths.get(".").resolve(path).toFile.getCanonicalPath)).toErrorOr flatMap { absolutePathToFile =>
      val file = File(absolutePathToFile)
      if (file.exists) {
        File(absolutePathToFile).contentAsString.validNel
      } else {
        s"Import file not found: $path".invalidNel
      }
    }
  }

  def specificFileResolver(filePath: Path): ImportResolver = CheckedAtoB.fromErrorOr { path =>
    if (path == filePath.toAbsolutePath.toString || path == filePath.name) {
      if (filePath.exists) {
        File(filePath.toAbsolutePath.pathAsString).contentAsString.validNel
      } else {
        s"Import file not found: $path".invalidNel
      }
    } else {
      s"'$path' did not match specific file ${filePath.toAbsolutePath.pathAsString}".invalidNel
    }
  }

  /**
    * Given a string representing an http(s) url, retrieve the contents
    */
  def httpResolver: ImportResolver =
    httpResolverWithHeaders(Map.empty[String,String])

  /**
    * Given a string representing an http(s) url, retrieve the contents
    * using the supplied headers (which can be an empty map)
    */
  def httpResolverWithHeaders(headers: Map[String,String]): ImportResolver = CheckedAtoB.fromCheck { str: String =>

    implicit val sttpBackend = AsyncHttpClientCatsBackend[IO]()

    val responseIO : IO[Response[String]] =
      sttp.get(uri"$str").headers(headers).send()

    // temporary situation to get functionality working before
    // starting in on async-ifying the entire WdlNamespace flow
    val result: Checked[String] = Await.result(responseIO.unsafeToFuture, 15.seconds).body.leftMap(NonEmptyList(_, List.empty))

    result
  }

}
