package cromwell.engine.backend.jes

import cromwell.engine.io.IoInterface
import wdl4s.values.WdlValue
import cromwell.engine._
import cromwell.util.PathUtil._
import cromwell.engine.Hashing._
import scala.util.Try

class JesWorkflowEngineFunctions(interface: IoInterface, context: WorkflowContext) extends WorkflowEngineFunctions(interface, context) {
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
