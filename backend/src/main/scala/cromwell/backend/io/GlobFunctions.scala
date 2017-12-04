package cromwell.backend.io

import cats.instances.list._
import cats.syntax.traverse._
import cromwell.backend.BackendJobDescriptor
import cromwell.core.CallContext
import common.validation.ErrorOr.ErrorOr
import wom.values._
import wom.expression.IoFunctionSet
import wom.graph.TaskCallNode
import wom.values.WomGlobFile

trait GlobFunctions extends IoFunctionSet {

  def callContext: CallContext

  def findGlobOutputs(call: TaskCallNode, jobDescriptor: BackendJobDescriptor): ErrorOr[List[WomGlobFile]] =
    call.callable.outputs.flatTraverse[ErrorOr, WomGlobFile] { outputDefinition =>
      outputDefinition.expression.evaluateFiles(jobDescriptor.localInputs, this, outputDefinition.womType) map {
        _.toList collect { case glob: WomGlobFile => glob }
      }
    }

  def globDirectory(glob: String): String = GlobFunctions.globName(glob) + "/"
  /**
    * Returns a path to the glob.
    *
    * This path is usually passed back into the glob() method below.
    *
    * @param glob The glob. This is the same "pattern" passed to glob() below.
    * @return The path.
    */
  def globPath(glob: String): String = callContext.root.resolve(globDirectory(glob)).pathAsString

  /**
    * Returns a list of path from the glob.
    *
    * The paths are read from a list file based on the pattern.
    *
    * @param pattern The pattern of the glob. This is the same "glob" passed to globPath().
    * @return The paths that match the pattern.
    */
  override def glob(pattern: String): Seq[String] = {
    import GlobFunctions._
    val globPatternName = globName(pattern)
    val listFilePath = callContext.root.resolve(s"${globName(pattern)}.list")
    // This "lines" is technically a read file and hence should use the readFile IO method
    listFilePath.toRealPath().lines.toList map { fileName =>
      (callContext.root /  globPatternName  / fileName).pathAsString
    }
  }
}

object GlobFunctions {
  def globName(glob: String) = s"glob-${glob.md5Sum}"
}
