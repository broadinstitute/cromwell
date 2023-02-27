package cromwell.backend.google.pipelines.batch


import cromwell.backend.google.pipelines.common.PipelinesApiExpressionFunctions
import cromwell.backend.standard.{StandardInitializationData, StandardValidatedRuntimeAttributesBuilder}

case class GcpBackendInitializationData(
                                         override val workflowPaths: GcpBatchWorkflowPaths,
                                         override val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder,
                                         gcpBatchConfiguration: GcpBatchConfiguration,
                                         //gcsCredentials: Credentials,
                                         //privateDockerEncryptionKeyName: Option[String],
                                         //privateDockerEncryptedToken: Option[String],

                                       ) extends StandardInitializationData(workflowPaths, runtimeAttributesBuilder, classOf[PipelinesApiExpressionFunctions] )
