package cromwell.engine.backend.jes

import cromwell.binding.values.WdlValue
import cromwell.engine._
import cromwell.engine.io.IoInterface
import cromwell.util.PathUtil._

import scala.util.Try
import Hashing._

class JesWorkflowEngineFunctions(interface: IoInterface, context: WorkflowContext) extends WorkflowEngineFunctions(interface, context){

  override def globPath(glob: String) = s"${context.root}/glob-${glob.md5Sum}/"

  override def fileContentsToString(path: String) = {
      if (!path.isGcsUrl) interface.readFile(s"${context.root}/$path")
      else interface.readFile(path)
  }
}

class JesCallEngineFunctions(interface: IoInterface, context: CallContext) extends JesWorkflowEngineFunctions(interface, context) with CallEngineFunctions {
  override def stdout(params: Seq[Try[WdlValue]]) = stdout(context)
  override def stderr(params: Seq[Try[WdlValue]]) = stderr(context)
}
