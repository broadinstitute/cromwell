package cromwell.backend.io

import cats.instances.list._
import cats.syntax.traverse._
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.BackendJobDescriptor
import cromwell.core.CallContext
import cromwell.core.io.AsyncIoFunctions
import wom.expression.IoFunctionSet
import wom.graph.CommandCallNode
import wom.values._

import scala.concurrent.Future

trait GlobFunctions extends IoFunctionSet with AsyncIoFunctions {

  def callContext: CallContext

  def findGlobOutputs(call: CommandCallNode, jobDescriptor: BackendJobDescriptor): ErrorOr[List[WomGlobFile]] = {
    def fromOutputs = call.callable.outputs.flatTraverse[ErrorOr, WomGlobFile] { outputDefinition =>
      outputDefinition.expression.evaluateFiles(jobDescriptor.localInputs, this, outputDefinition.womType) map {
        _.toList.flatMap(_.file.flattenFiles) collect { case glob: WomGlobFile => glob }
      }
    }
    fromOutputs.map(_ ++ call.callable.additionalGlob)
  }

  /**
    * Returns a list of path from the glob.
    *
    * The paths are read from a list file based on the pattern.
    *
    * @param pattern The pattern of the glob. This is the same "glob" passed to globPath().
    * @return The paths that match the pattern.
    */
  override def glob(pattern: String): Future[Seq[String]] = {
    import GlobFunctions._
    val globPatternName = globName(pattern)
    val listFilePath = callContext.root.resolve(s"${globName(pattern)}.list")
    asyncIo.readLinesAsync(listFilePath.toRealPath()) map { lines =>
      lines.toList map { fileName =>
        (callContext.root /  globPatternName  / fileName).pathAsString
      }
    }
  }
}
