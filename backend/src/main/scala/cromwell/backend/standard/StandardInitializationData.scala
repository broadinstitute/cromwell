package cromwell.backend.standard

import cromwell.backend.BackendInitializationData
import cromwell.backend.io.WorkflowPaths
import cromwell.core.path.PathBuilder

class StandardInitializationData
(
  val workflowPaths: WorkflowPaths,
  val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder
) extends BackendInitializationData

object StandardInitializationData {
  def workflowPaths(initializationDataOption: Option[BackendInitializationData]): WorkflowPaths = {
    BackendInitializationData.as[StandardInitializationData](initializationDataOption).workflowPaths
  }

  def pathBuilders(initializationDataOption: Option[BackendInitializationData]): List[PathBuilder] = {
    workflowPaths(initializationDataOption).pathBuilders
  }
}
