package cromwell.backend.io

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.backend.io.DirectoryFunctions._
import cromwell.core.path.{Path, PathFactory}
import wom.expression.IoFunctionSet
import wom.values.{WomFile, WomGlobFile, WomMaybeListedDirectory, WomMaybePopulatedFile, WomSingleFile, WomUnlistedDirectory}

import scala.concurrent.Future

trait DirectoryFunctions extends IoFunctionSet with PathFactory {
  override def listAllFilesUnderDirectory(dirPath: String): Future[Seq[String]] = {
    temporaryImplListPaths(dirPath)
  }

  // TODO: WOM: WOMFILE: This will likely use a Tuple2(tar file, dir list file) for each dirPath.
  private final def temporaryImplListPaths(dirPath: String): Future[Seq[String]] = {
    val errorOrPaths = for {
      dir <- validate(buildPath(ensureSlashed(dirPath)))
      files <- listFiles(dir)
    } yield files.map(_.pathAsString)
    Future.fromTry(errorOrPaths.toTry(s"Error listing files under $dirPath"))
  }
}

object DirectoryFunctions {
  def ensureSlashed(dir: String): String = if (dir.endsWith("/")) dir else s"$dir/"

  def listFiles(path: Path): ErrorOr[List[Path]] = {
    if (path.isDirectory) {
      for {
        pathListing <- validate(path.list.toList)
        pathsPerListing <- pathListing.traverse(listFiles)
      } yield pathsPerListing.flatten
    } else {
      List(path).valid
    }
  }

  def listWomSingleFiles(womFile: WomFile, pathFactory: PathFactory, pathPatcher: String => String): ErrorOr[List[WomSingleFile]] = {
    def listWomSingleFiles(womFile: WomFile): ErrorOr[List[WomSingleFile]] = {
      womFile match {
        case womSingleFile: WomSingleFile => List(womSingleFile).valid

        case womUnlistedDirectory: WomUnlistedDirectory =>
          val errorOrListPaths = listFiles(pathFactory.buildPath(ensureSlashed(pathPatcher(womUnlistedDirectory.value))))
          errorOrListPaths.map(_.map(path => WomSingleFile(path.pathAsString)))

        case womMaybePopulatedFile: WomMaybePopulatedFile =>
          val allFiles: List[WomFile] =
            womMaybePopulatedFile.valueOption.toList.map(WomSingleFile) ++ womMaybePopulatedFile.secondaryFiles
          allFiles.traverse(listWomSingleFiles).map(_.flatten)

        case w: WomMaybeListedDirectory =>
          (w.valueOption, w.listingOption) match {
            case (None, None) => Nil.valid
            case (Some(value), None) => listWomSingleFiles(WomUnlistedDirectory(value))
            case (None, Some(listing)) => listing.toList.traverse(listWomSingleFiles).map(_.flatten)
            // TODO: WOM: WOMFILE: This is a special case where files from a different path are supposed to end up under the directory. If this implementation is correct, remove this TODO.
            case (Some(_), Some(listing)) => listing.toList.traverse(listWomSingleFiles).map(_.flatten)
          }
        // TODO: WOM: WOMFILE: How did a glob get here? Should this link into glob functions to list the globs?

        case _: WomGlobFile => s"Unexpected glob / unable to list glob files at this time: $womFile".invalidNel
      }
    }

    listWomSingleFiles(womFile)
  }
}
