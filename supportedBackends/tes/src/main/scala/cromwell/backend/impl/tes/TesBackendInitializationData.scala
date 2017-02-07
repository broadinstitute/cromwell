package cromwell.backend.impl.tes

import cromwell.backend.standard.{StandardExpressionFunctions, StandardInitializationData, StandardValidatedRuntimeAttributesBuilder}

case class TesBackendInitializationData
(
  override val workflowPaths: TesWorkflowPaths,
  override val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder,
  tesConfiguration: TesConfiguration
) extends StandardInitializationData(workflowPaths, runtimeAttributesBuilder, classOf[StandardExpressionFunctions])
