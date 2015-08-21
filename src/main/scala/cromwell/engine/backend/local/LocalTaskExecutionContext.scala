package cromwell.engine.backend.local

import java.nio.file.{Path, Paths}

import cromwell.engine.backend.TaskExecutionContext

case class LocalTaskExecutionContext(workflowRoot: Path, stdout: Path, stderr: Path, cwd: Path = Paths.get(".")) extends TaskExecutionContext {
  override def engineFunctions: LocalEngineFunctions = new LocalEngineFunctions(this)
}
