package cromwell.engine.backend.local

import java.nio.file.Paths
import wdl4s.values.WdlValue
import cromwell.engine.io.IoInterface
import cromwell.engine.{CallContext, CallEngineFunctions, WorkflowContext, WorkflowEngineFunctions, _}

import scala.language.postfixOps
import scala.util.Try

class LocalWorkflowEngineFunctions(interface: IoInterface, context: WorkflowContext) extends WorkflowEngineFunctions(interface, context)

class LocalCallEngineFunctions(interface: IoInterface, context: CallContext) extends LocalWorkflowEngineFunctions(interface ,context) with CallEngineFunctions {
  import PathString._

  override def adjustFilePath(path: String) = {
    if (!Paths.get(path).isAbsolute && !path.isUriWithProtocol) Paths.get(context.root).resolve(path).toString else path
  }

  override def stdout(params: Seq[Try[WdlValue]]) = stdout(context)
  override def stderr(params: Seq[Try[WdlValue]]) = stderr(context)
}

