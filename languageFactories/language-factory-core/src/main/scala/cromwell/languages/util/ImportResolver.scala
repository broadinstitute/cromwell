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
import common.validation.Checked._
import common.validation.Validation._
import cromwell.core.path.Path
import wom.core.WorkflowSource

import scala.concurrent.duration._
import scala.concurrent.Await
import scala.util.{Failure, Success, Try}

object ImportResolver {
  trait ImportResolver2 {
    def name: String
    protected def innerResolver(str: String): Checked[WorkflowSource]
    def updateForNewRoot(newRoot: String): ImportResolver2
    def resolver: CheckedAtoB[String, WorkflowSource] = CheckedAtoB.fromCheck { str =>
      innerResolver(str).contextualizeErrors(s"resolve '$str' using resolver: '$name'")
    }
  }

  object DirectoryResolver2 {
    def apply(directory: Path, allowEscapingDirectory: Boolean): DirectoryResolver2 = {
      val dontEscapeFrom = if (allowEscapingDirectory) None else Option(directory.toJava.getCanonicalPath)
      DirectoryResolver2(directory, dontEscapeFrom)
    }
  }

  case class DirectoryResolver2(directory: Path, dontEscapeFrom: Option[String] = None) extends ImportResolver2 {
    lazy val absolutePathToDirectory = directory.toJava.getCanonicalPath

    override def innerResolver(path: String): Checked[WorkflowSource] = {
      Try(Paths.get(directory.resolve(path).toFile.getCanonicalPath)).toErrorOr flatMap { absolutePathToFile =>
        if (dontEscapeFrom.forall(absolutePathToFile.startsWith)) {
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
    }.toEither

    override def updateForNewRoot(newRoot: String) = {
      val newParentDirectory = directory.resolve(newRoot).parent
      DirectoryResolver2(newParentDirectory, dontEscapeFrom)
    }

    override def name: String = s"relative to directory $absolutePathToDirectory (without escaping $dontEscapeFrom)"
  }

  def zippedImportResolver2(zippedImports: Array[Byte]): ErrorOr[ImportResolver2] = {
    LanguageFactoryUtil.validateImportsDirectory(zippedImports) map { dir =>
      DirectoryResolver2(dir, Option(dir.toJava.getCanonicalPath))
    }
  }

  case object AnyLocalFileResolver2 extends ImportResolver2 {
    override def innerResolver(path: String): Checked[WorkflowSource] = {
      Try(Paths.get(Paths.get(".").resolve(path).toFile.getCanonicalPath)).toErrorOr flatMap { absolutePathToFile =>
        val file = File(absolutePathToFile)
        if (file.exists) {
          File(absolutePathToFile).contentAsString.validNel
        } else {
          s"Import file not found: $path".invalidNel
        }
      }
    }.toEither
    override def updateForNewRoot(newRoot: String): ImportResolver2 = this
    override val name = "any local file"
  }

  case class SpecificFileResolver2(filePath: Path) extends ImportResolver2 {
    override def innerResolver(path: String): Checked[WorkflowSource] = {
      if (path == filePath.toAbsolutePath.toString || path == filePath.name) {
        if (filePath.exists) {
          File(filePath.toAbsolutePath.pathAsString).contentAsString.validNelCheck
        } else {
          s"Import file not found: $path".invalidNelCheck
        }
      } else {
        s"'$path' did not match".invalidNelCheck
      }
    }

    override def updateForNewRoot(newRoot: String): ImportResolver2 = this
    override lazy val name = s"specific file ${filePath.toJava.getCanonicalPath}"
  }

  case class HttpResolver2(relativeTo: Option[String] = None, headers: Map[String, String] = Map.empty) extends ImportResolver2 {
    override def name: String = {
      s"http importer${relativeTo map { r => s" (relative to $r)" } getOrElse ""}"
    }

    override def innerResolver(str: String): Checked[WorkflowSource] = {
      val toLookup = if (str.startsWith("http")) str else s"$relativeTo/$str"

      implicit val sttpBackend = AsyncHttpClientCatsBackend[IO]()
      Try {
        val responseIO: IO[Response[String]] = sttp.get(uri"$toLookup").headers(headers).send()

        // temporary situation to get functionality working before
        // starting in on async-ifying the entire WdlNamespace flow
        val result: Checked[String] = Await.result(responseIO.unsafeToFuture, 15.seconds).body.leftMap(NonEmptyList(_, List.empty))

        result
      } match {
        case Success(result) => result
        case Failure(e) => s"HTTP resolver with headers had an unexpected error (${e.getMessage})".invalidNelCheck
      }
    }

    override def updateForNewRoot(newRoot: String): ImportResolver2 = HttpResolver2(Option(newRoot), headers)
  }

}
