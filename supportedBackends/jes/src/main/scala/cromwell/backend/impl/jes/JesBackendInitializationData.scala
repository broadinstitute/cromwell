package cromwell.backend.impl.jes

import com.google.api.services.genomics.Genomics
import cromwell.backend.standard.{StandardInitializationData, StandardValidatedRuntimeAttributesBuilder}

case class JesBackendInitializationData
(
  override val workflowPaths: JesWorkflowPaths,
  override val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder,
  jesConfiguration: JesConfiguration,
  genomics: Genomics
) extends StandardInitializationData(workflowPaths, runtimeAttributesBuilder, classOf[JesExpressionFunctions])
