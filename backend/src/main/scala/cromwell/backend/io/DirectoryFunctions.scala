package cromwell.backend.io

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.util.StringUtil._
import common.validation.ErrorOr._
import cromwell.backend.BackendJobDescriptor
import cromwell.core.path.{Path, PathFactory}
import wom.expression.IoFunctionSet
import wom.expression.IoFunctionSet.{IoDirectory, IoElement, IoFile}
import wom.graph.CommandCallNode
import wom.values.{WomFile, WomGlobFile, WomMaybeListedDirectory, WomMaybePopulatedFile, WomSingleFile, WomUnlistedDirectory}

import scala.concurrent.Future
import scala.util.Try

trait DirectoryFunctions extends IoFunctionSet with PathFactory {

  def findDirectoryOutputs(call: CommandCallNode,
                           jobDescriptor: BackendJobDescriptor): ErrorOr[List[WomUnlistedDirectory]] = {
    call.callable.outputs.flatTraverse { outputDefinition =>
      outputDefinition.expression.evaluateFiles(jobDescriptor.localInputs, this, outputDefinition.womType) map {
        _.toList.flatMap(_.flattenFiles) collect { case unlistedDirectory: WomUnlistedDirectory => unlistedDirectory }
      }
    }
  }

  override def isDirectory(path: String) = Future.fromTry(Try(buildPath(path).isDirectory))

  override def listDirectory(path: String)(visited: Vector[String] = Vector.empty): Future[Iterator[IoElement]] = {
    Future.fromTry(Try {
      val visitedPaths = visited.map(buildPath)
      val cromwellPath = buildPath(path.ensureSlashed)

      // To prevent infinite recursion through symbolic links make sure we don't visit the same directory twice
      def hasBeenVisited(other: Path) = visitedPaths.exists(_.isSameFileAs(other))

      cromwellPath.list.collect({
        case directory if directory.isDirectory &&
          !cromwellPath.isSamePathAs(directory) &&
          !hasBeenVisited(directory) => IoDirectory(directory.pathAsString)
        case file => IoFile(file.pathAsString)
      })
    })
  }
}

object DirectoryFunctions {
  def listFiles(path: Path): ErrorOr[List[Path]] = path.listRecursively.filterNot(_.isDirectory).toList.validNel

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
