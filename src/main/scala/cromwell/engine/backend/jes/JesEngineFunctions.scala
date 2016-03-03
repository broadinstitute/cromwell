package cromwell.engine.backend.jes

import cromwell.engine.Hashing._
import cromwell.engine.PathString._
import cromwell.engine._
import cromwell.engine.io.IoInterface
import wdl4s.values.WdlValue

import scala.util.Try
class JesWorkflowEngineFunctions(interface: IoInterface, context: WorkflowContext) extends WorkflowEngineFunctions(interface, context) {
  override def globPath(glob: String) = s"${context.root}/glob-${glob.md5Sum}/"

  override def adjustFilePath(path: String) = if (!path.isGcsUrl) s"${context.root}/$path" else path
}

class JesCallEngineFunctions(interface: IoInterface, context: CallContext) extends JesWorkflowEngineFunctions(interface, context) with CallEngineFunctions {
  override def stdout(params: Seq[Try[WdlValue]]) = stdout(context)
  override def stderr(params: Seq[Try[WdlValue]]) = stderr(context)
}
