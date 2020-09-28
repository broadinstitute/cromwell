package cromwell.backend.google.pipelines.common

import com.google.auth.Credentials
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory
import cromwell.backend.standard.{StandardInitializationData, StandardValidatedRuntimeAttributesBuilder}

case class PipelinesApiBackendInitializationData
(
  override val workflowPaths: PipelinesApiWorkflowPaths,
  override val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder,
  papiConfiguration: PipelinesApiConfiguration,
  gcsCredentials: Credentials,
  genomicsRequestFactory: PipelinesApiRequestFactory,
  privateDockerEncryptionKeyName: Option[String],
  privateDockerEncryptedToken: Option[String],
  vpcNetworkAndSubnetworkProjectLabels: Option[VpcAndSubnetworkProjectLabelValues]
) extends StandardInitializationData(workflowPaths, runtimeAttributesBuilder, classOf[PipelinesApiExpressionFunctions])
