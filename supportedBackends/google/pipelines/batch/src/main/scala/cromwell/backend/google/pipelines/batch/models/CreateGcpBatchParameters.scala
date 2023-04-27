package cromwell.backend.google.pipelines.batch.models

import cromwell.backend.BackendJobDescriptor

case class CreateGcpBatchParameters(jobDescriptor: BackendJobDescriptor,
                                    runtimeAttributes: GcpBatchRuntimeAttributes,
                                    batchAttributes: GcpBatchConfigurationAttributes,
                                    projectId: String,
                                    region: String
                                   )
