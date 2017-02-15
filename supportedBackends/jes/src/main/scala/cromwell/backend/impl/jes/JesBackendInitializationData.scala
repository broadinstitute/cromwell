package cromwell.backend.impl.jes

import com.google.api.services.genomics.Genomics
import com.google.auth.Credentials
import cromwell.backend.standard.{StandardInitializationData, StandardValidatedRuntimeAttributesBuilder}

case class JesBackendInitializationData
(
  override val workflowPaths: JesWorkflowPaths,
  override val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder,
  jesConfiguration: JesConfiguration,
  gcsCredentials: Credentials,
  genomics: Genomics
) extends StandardInitializationData(workflowPaths, runtimeAttributesBuilder, classOf[JesExpressionFunctions])
