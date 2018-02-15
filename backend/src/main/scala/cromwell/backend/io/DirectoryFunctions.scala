package cromwell.backend.io

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.util.StringUtil._
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.io.DirectoryFunctions._
import cromwell.core.path.{Path, PathFactory}
import wom.expression.IoFunctionSet
import wom.graph.CommandCallNode
import wom.values.{WomFile, WomGlobFile, WomMaybeListedDirectory, WomMaybePopulatedFile, WomSingleFile, WomUnlistedDirectory}

import scala.concurrent.Future

trait DirectoryFunctions extends IoFunctionSet with PathFactory {

  def findDirectoryOutputs(call: CommandCallNode,
                           jobDescriptor: BackendJobDescriptor): ErrorOr[List[WomUnlistedDirectory]] = {
    call.callable.outputs.flatTraverse[ErrorOr, WomUnlistedDirectory] { outputDefinition =>
      outputDefinition.expression.evaluateFiles(jobDescriptor.localInputs, this, outputDefinition.womType) map {
        _.toList.flatMap(_.flattenFiles) collect { case unlistedDirectory: WomUnlistedDirectory => unlistedDirectory }
      }
    }
  }

  override def listAllFilesUnderDirectory(dirPath: String): Future[Seq[String]] = {
    temporaryImplListPaths(dirPath)
  }

  // TODO: WOM: WOMFILE: This will likely use a Tuple2(tar file, dir list file) for each dirPath.
  private final def temporaryImplListPaths(dirPath: String): Future[Seq[String]] = {
    val errorOrPaths = for {
      dir <- validate(buildPath(dirPath.ensureSlashed))
      files <- listFiles(dir)
    } yield files.map(_.pathAsString)
    Future.fromTry(errorOrPaths.toTry(s"Error listing files under $dirPath"))
  }
}

object DirectoryFunctions {
  def listFiles(path: Path): ErrorOr[List[Path]] = {
    def listPaths(path: Path, checkedPaths: Set[Path]): ErrorOr[Set[Path]] = {
      val newCheckedPaths = checkedPaths ++ Set(path)
      if (path.isDirectory) {
        for {
          pathListing <- validate(path.list.toSet)
          uncheckedPaths = pathListing -- newCheckedPaths
          pathsPerListing <-
            uncheckedPaths.toList.traverse[ErrorOr, List[Path]](listPaths(_, newCheckedPaths).map(_.toList))
        } yield checkedPaths ++ pathsPerListing.flatten.toSet
      } else {
        newCheckedPaths.valid
      }
    }

    val allPaths = listPaths(path, Set.empty).map(_.toList)
    val allFiles = allPaths.map(_.filterNot(_.isDirectory))
    allFiles
  }

  def listWomSingleFiles(womFile: WomFile, pathFactory: PathFactory): ErrorOr[List[WomSingleFile]] = {
    def listWomSingleFiles(womFile: WomFile): ErrorOr[List[WomSingleFile]] = {
      womFile match {
        case womSingleFile: WomSingleFile => List(womSingleFile).valid

        case womUnlistedDirectory: WomUnlistedDirectory =>
          val errorOrListPaths = listFiles(pathFactory.buildPath(womUnlistedDirectory.value.ensureSlashed))
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
