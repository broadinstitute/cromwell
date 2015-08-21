package cromwell.engine.backend.jes

import cromwell.engine.backend.TaskExecutionContext
import cromwell.util.google.GoogleCloudStoragePath

case class JesTaskExecutionContext(callDir: GoogleCloudStoragePath, jesConnection: JesInterface) extends TaskExecutionContext {
  override def engineFunctions: JesEngineFunctions = new JesEngineFunctions(this)
}
