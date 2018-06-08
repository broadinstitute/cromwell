package cromwell.backend.io

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.util.StringUtil._
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.io.DirectoryFunctions.listFiles
import cromwell.core.io.AsyncIoFunctions
import cromwell.core.path.{Path, PathFactory}
import wom.expression.IoFunctionSet
import wom.expression.IoFunctionSet.{IoDirectory, IoElement, IoFile}
import wom.graph.CommandCallNode
import wom.values.{WomFile, WomGlobFile, WomMaybeListedDirectory, WomMaybePopulatedFile, WomSingleFile, WomUnlistedDirectory}

import scala.concurrent.Future
import scala.util.Try

trait DirectoryFunctions extends IoFunctionSet with PathFactory with AsyncIoFunctions {

  def findDirectoryOutputs(call: CommandCallNode,
                           jobDescriptor: BackendJobDescriptor): ErrorOr[List[WomUnlistedDirectory]] = {
    call.callable.outputs.flatTraverse[ErrorOr, WomUnlistedDirectory] { outputDefinition =>
      outputDefinition.expression.evaluateFiles(jobDescriptor.localInputs, this, outputDefinition.womType) map {
        _.toList.flatMap(_.file.flattenFiles) collect { case unlistedDirectory: WomUnlistedDirectory => unlistedDirectory }
      }
    }
  }

  override def isDirectory(path: String) = asyncIo.isDirectory(buildPath(path))

  /*
   * Several things are wrong here.
   * 1) None of this is going through the I/O Actor: https://github.com/broadinstitute/cromwell/issues/3133
   * which means no instrumentation, no throttling, no batching, and no custom retries.
   * 2) The NIO implementation of "list" in GCS will list all objects with the prefix "path", unlike the unix
   * implementation which lists files and directories children. What we need is the unix behavior, even for cloud filesystems.
   * 3) It uses the isDirectory function directly on the path, which cannot be trusted for GCS paths. It should use asyncIo.isDirectory instead.
   */
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
