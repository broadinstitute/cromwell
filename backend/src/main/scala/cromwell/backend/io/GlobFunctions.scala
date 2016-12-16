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

  override def globPath(glob: String): String = callContext.root.resolve(globDirectory(glob)).toString

  override def glob(path: String, pattern: String): Seq[String] = {
    val name = globName(pattern)
    val listFile = callContext.root.resolve(s"$name.list").toRealPath()
    Files.readAllLines(listFile).asScala map { fileName => callContext.root.resolve(s"$name/$fileName").toRealString }
  }
}
