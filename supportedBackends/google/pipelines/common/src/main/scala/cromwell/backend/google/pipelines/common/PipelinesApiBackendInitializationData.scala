package cromwell.backend.google.pipelines.common

import com.google.auth.Credentials
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory
import cromwell.backend.standard.{StandardInitializationData, StandardValidatedRuntimeAttributesBuilder}

case class PipelinesApiBackendInitializationData
(
  override val workflowPaths: PipelinesApiWorkflowPaths,
  override val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder,
  jesConfiguration: PipelinesApiConfiguration,
  gcsCredentials: Credentials,
  genomicsRequestFactory: PipelinesApiRequestFactory
) extends StandardInitializationData(workflowPaths, runtimeAttributesBuilder, classOf[PipelinesApiExpressionFunctions])
