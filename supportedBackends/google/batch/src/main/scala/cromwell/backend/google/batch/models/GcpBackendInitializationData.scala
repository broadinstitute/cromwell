package cromwell.backend.google.batch.models

import com.google.auth.Credentials
import cromwell.backend.google.batch.api.GcpBatchRequestFactory
import cromwell.backend.google.batch.util.BatchExpressionFunctions
import cromwell.backend.standard.{StandardInitializationData, StandardValidatedRuntimeAttributesBuilder}

case class GcpBackendInitializationData(
  override val workflowPaths: GcpBatchWorkflowPaths,
  override val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder,
  gcpBatchConfiguration: GcpBatchConfiguration,
  gcsCredentials: Credentials,
  privateDockerEncryptionKeyName: Option[String],
  privateDockerEncryptedToken: Option[String],
  vpcNetworkAndSubnetworkProjectLabels: Option[VpcAndSubnetworkProjectLabelValues],
  requestFactory: GcpBatchRequestFactory
) extends StandardInitializationData(workflowPaths, runtimeAttributesBuilder, classOf[BatchExpressionFunctions])
