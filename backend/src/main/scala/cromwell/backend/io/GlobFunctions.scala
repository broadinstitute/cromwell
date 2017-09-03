package cromwell.backend.io

import cromwell.backend.BackendJobDescriptor
import cromwell.core.CallContext
import wdl4s.wdl.values._
import wdl4s.wom.expression.IoFunctionSet
import wdl4s.wom.graph.TaskCallNode

trait GlobFunctions extends IoFunctionSet {

  def callContext: CallContext

  def findGlobOutputs(call: TaskCallNode, jobDescriptor: BackendJobDescriptor): Set[WdlGlobFile] = {
    // TODO WOM: How to get all the WdlGlobFile ?
    val globOutputs = Set.empty[WdlGlobFile]

    globOutputs
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
    val name = globName(pattern)
    val listFile = callContext.root.resolve(s"$name.list").toRealPath()
    listFile.lines.toSeq map { fileName =>
      callContext.root.resolve(name).resolve(fileName).pathAsString
    }
  }
}
