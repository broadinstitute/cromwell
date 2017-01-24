package cromwell.backend.impl.tes

import cromwell.backend.io.WorkflowPathsWithDocker
import cromwell.backend.standard.{StandardInitializationData, StandardValidatedRuntimeAttributesBuilder}

case class TesBackendInitializationData
(
  override val workflowPaths: WorkflowPathsWithDocker,
  override val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder,
  tesConfiguration: TesConfiguration
) extends StandardInitializationData(workflowPaths, runtimeAttributesBuilder)