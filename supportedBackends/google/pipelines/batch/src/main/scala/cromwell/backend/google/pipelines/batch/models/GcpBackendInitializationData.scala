package cromwell.backend.google.pipelines.batch.models

import com.google.auth.Credentials
import cromwell.backend.google.pipelines.batch.util.BatchExpressionFunctions
import cromwell.backend.standard.{StandardInitializationData, StandardValidatedRuntimeAttributesBuilder}

case class GcpBackendInitializationData(
                                         override val workflowPaths: GcpBatchWorkflowPaths,
                                         override val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder,
                                         gcpBatchConfiguration: GcpBatchConfiguration,
                                         gcsCredentials: Credentials,
                                         privateDockerEncryptionKeyName: Option[String],
                                         privateDockerEncryptedToken: Option[String],
                                         vpcNetworkAndSubnetworkProjectLabels: Option[VpcAndSubnetworkProjectLabelValues]

                                       ) extends StandardInitializationData(workflowPaths, runtimeAttributesBuilder, classOf[BatchExpressionFunctions] )
