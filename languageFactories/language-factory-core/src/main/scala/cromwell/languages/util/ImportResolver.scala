package cromwell.languages.util

import java.nio.file.Paths
import better.files.File
import cats.syntax.either._
import cats.syntax.validated._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.core.path.Path
import wom.core.{WorkflowSource}
import scala.util.Try

object ImportResolver {
  type ImportResolver = CheckedAtoB[String, WorkflowSource]

  def zippedImportsResolver(zippedImports: Array[Byte]): ImportResolver = {
    directoryResolver(LanguageFactoryUtil.validateImportsDirectory(zippedImports).toEither)
  }

  def directoryResolver(directory: Path): ImportResolver = CheckedAtoB.fromErrorOr { path =>
    Try(Paths.get(directory.resolve(path).toFile.getCanonicalPath)).toErrorOr flatMap { absolutePathToFile =>
      val absolutePathToImports = Paths.get(directory.toJava.getCanonicalPath)
      if (absolutePathToFile.startsWith(absolutePathToImports)) {
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

}
