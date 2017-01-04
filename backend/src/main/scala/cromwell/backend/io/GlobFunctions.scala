package cromwell.backend.io

import java.nio.file.Files

import cromwell.backend.BackendJobDescriptor
import cromwell.core.CallContext
import wdl4s.TaskCall
import wdl4s.expression.{NoFunctions, PureStandardLibraryFunctionsLike}
import wdl4s.values._
import scala.collection.JavaConverters._
import cromwell.core.path.PathImplicits._

trait GlobFunctions extends PureStandardLibraryFunctionsLike {

  def callContext: CallContext

  def findGlobOutputs(call: TaskCall, jobDescriptor: BackendJobDescriptor): Set[WdlGlobFile] = {
    val globOutputs = call.task.findOutputFiles(jobDescriptor.fullyQualifiedInputs, NoFunctions) collect {
      case glob: WdlGlobFile => glob
    }

    globOutputs.distinct.toSet
  }

  def globDirectory(glob: String): String = globName(glob) + "/"
  def globName(glob: String) = s"glob-${glob.md5Sum}"

  /**
    * Returns a path to the glob using toRealString.
    *
    * NOTE: Due to use of toRealString, returned paths must be passed to PathBuilders.buildPath, and will not work with
    * Paths.get.
    *
    * This path is usually passed back into the glob() method below.
    *
    * @param glob The glob. This is the same "pattern" passed to glob() below.
    * @return The path converted using .toRealString.
    */
  override def globPath(glob: String): String = callContext.root.resolve(globDirectory(glob)).toRealString

  /**
    * Returns a list of path from the glob, each path converted to a string using toRealString.
    *
    * The paths are currently read from a list file based on the pattern, and the path parameter is not used.
    *
    * NOTE: Due to use of toRealString, returned paths must be passed to PathBuilders.buildPath, and will not work with
    * Paths.get.
    *
    * @param path    The path string returned by globPath. This isn't currently used.
    * @param pattern The pattern of the glob. This is the same "glob" passed to globPath().
    * @return The paths that match the pattern, each path converted using .toRealString.
    */
  override def glob(path: String, pattern: String): Seq[String] = {
    val name = globName(pattern)
    val listFile = callContext.root.resolve(s"$name.list").toRealPath()
    Files.readAllLines(listFile).asScala map { fileName =>
      callContext.root.resolve(name).resolve(fileName).toRealString
    }
  }
}
