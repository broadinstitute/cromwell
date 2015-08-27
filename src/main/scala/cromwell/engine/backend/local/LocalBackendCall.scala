package cromwell.engine.backend.local

import java.nio.file.Paths

import cromwell.binding.values.WdlValue
import cromwell.binding.{CallOutputs, CallInputs, Call, WorkflowDescriptor}
import cromwell.engine.AbortRegistrationFunction
import cromwell.engine.backend.{LocalFileSystemBackendCall, BackendCall}

import scala.util.Try

case class LocalBackendCall(backend: LocalBackend,
                            workflowDescriptor: WorkflowDescriptor,
                            call: Call,
                            locallyQualifiedInputs: CallInputs,
                            callAbortRegistrationFunction: AbortRegistrationFunction) extends BackendCall with LocalFileSystemBackendCall {
  val workflowRootPath = LocalBackend.hostExecutionPath(workflowDescriptor)
  val callRootPath = LocalBackend.hostCallPath(workflowDescriptor, call.name)
  val dockerContainerExecutionDir = s"/root/${workflowDescriptor.id.toString}"
  val containerCallRoot = call.docker match {
    case Some(docker) => Paths.get(dockerContainerExecutionDir).resolve(s"call-${call.name}")
    case None => callRootPath
  }
  val rc = callRootPath.resolve("rc")
  val stdout = callRootPath.resolve("stdout")
  val stderr = callRootPath.resolve("stderr")
  val script = callRootPath.resolve("script")
  val engineFunctions: LocalEngineFunctions = new LocalEngineFunctions(callRootPath, stdout, stderr)
  val lookupFunction: String => WdlValue = inputName => locallyQualifiedInputs.get(inputName).get
  callRootPath.toFile.mkdirs
  def execute: Try[CallOutputs] = backend.execute(this)
}
