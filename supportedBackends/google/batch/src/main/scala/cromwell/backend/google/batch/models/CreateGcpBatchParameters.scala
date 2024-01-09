package cromwell.backend.google.batch.models

import cromwell.backend.BackendJobDescriptor

case class CreateGcpBatchParameters(jobDescriptor: BackendJobDescriptor,
                                    runtimeAttributes: GcpBatchRuntimeAttributes,
                                    batchAttributes: GcpBatchConfigurationAttributes,
                                    projectId: String,
                                    region: String
)
