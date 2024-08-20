package cromwell.backend.io

import cats.implicits._
import common.util.StringUtil._
import common.validation.ErrorOr._
import cromwell.backend.BackendJobDescriptor
import cromwell.core.io.AsyncIoFunctions
import cromwell.core.path.{Path, PathFactory}
import wom.expression.{IoFunctionSet, IoFunctionSetAdapter}
import wom.graph.CommandCallNode
import wom.values.{WomFile, WomGlobFile, WomSingleFile, WomUnlistedDirectory}

trait DirectoryFunctions extends IoFunctionSet with PathFactory with AsyncIoFunctions {

  private lazy val evaluateFileFunctions = new IoFunctionSetAdapter(this) with FileEvaluationIoFunctionSet

  def findDirectoryOutputs(call: CommandCallNode,
                           jobDescriptor: BackendJobDescriptor
  ): ErrorOr[List[WomUnlistedDirectory]] =
    call.callable.outputs.flatTraverse[ErrorOr, WomUnlistedDirectory] { outputDefinition =>
      outputDefinition.expression.evaluateFiles(jobDescriptor.localInputs,
                                                evaluateFileFunctions,
                                                outputDefinition.womType
      ) map {
        _.toList.flatMap(_.file.flattenFiles) collect { case unlistedDirectory: WomUnlistedDirectory =>
          unlistedDirectory
        }
      }
    }
}

object DirectoryFunctions {
  private def listFiles(path: Path): ErrorOr[List[Path]] = path.listRecursively.filterNot(_.isDirectory).toList.validNel

  def listWomSingleFiles(womFile: WomFile, pathFactory: PathFactory): ErrorOr[List[WomSingleFile]] = {
    def listWomSingleFiles(womFile: WomFile): ErrorOr[List[WomSingleFile]] =
      womFile match {
        case womSingleFile: WomSingleFile => List(womSingleFile).valid

        case womUnlistedDirectory: WomUnlistedDirectory =>
          val errorOrListPaths = listFiles(pathFactory.buildPath(womUnlistedDirectory.value.ensureSlashed))
          errorOrListPaths.map(_.map(path => WomSingleFile(path.pathAsString)))

        // TODO: WOM: WOMFILE: How did a glob get here? Should this link into glob functions to list the globs?
        case _: WomGlobFile => s"Unexpected glob / unable to list glob files at this time: $womFile".invalidNel
      }

    listWomSingleFiles(womFile)
  }
}
