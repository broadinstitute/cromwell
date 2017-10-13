package cromwell.backend.io

import cats.instances.list._
import cats.syntax.traverse._
import cromwell.backend.BackendJobDescriptor
import cromwell.core.CallContext
import lenthall.validation.ErrorOr.ErrorOr
import wom.values._
import wom.expression.IoFunctionSet
import wom.graph.TaskCallNode
import wom.types.WdlAnyType
import wom.values.WdlGlobFile

trait GlobFunctions extends IoFunctionSet {

  def callContext: CallContext

  def findGlobOutputs(call: TaskCallNode, jobDescriptor: BackendJobDescriptor): ErrorOr[List[WdlGlobFile]] =
    call.callable.outputs.flatTraverse[ErrorOr, WdlGlobFile] {
      _.expression.evaluateFiles(jobDescriptor.fullyQualifiedInputs, this, WdlAnyType) map {
        _.toList collect { case glob: WdlGlobFile => glob }
      }
    }

  def globDirectory(glob: String): String = globName(glob) + "/"
  def globName(glob: String) = s"glob-${glob.md5Sum}"

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
    * The paths are currently read from a list file based on the pattern, and the path parameter is not used.
    *
    * @param path    The path string returned by globPath. This isn't currently used.
    * @param pattern The pattern of the glob. This is the same "glob" passed to globPath().
    * @return The paths that match the pattern.
    */
  override def glob(path: String, pattern: String): Seq[String] = {
    val globPatternName = globName(pattern)
    val listFilePath = callContext.root.resolve(s"${globName(pattern)}.list")
    // This "lines" is technically a read file and hence should use the readFile IO method
    listFilePath.toRealPath().lines.toList map { fileName =>
      (callContext.root /  globPatternName  / fileName).pathAsString
    }
  }
}
