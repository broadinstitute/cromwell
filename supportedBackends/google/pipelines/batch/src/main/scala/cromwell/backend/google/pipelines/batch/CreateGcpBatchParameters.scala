package cromwell.backend.google.pipelines.batch

import cromwell.backend.BackendJobDescriptor

case class CreateGcpBatchParameters(jobDescriptor: BackendJobDescriptor,
                                    runtimeAttributes: GcpBatchRuntimeAttributes,
                                    batchAttributes: GcpBatchConfigurationAttributes,
                                    dockerImage: String,
                                    projectId: String,
                                    region: String
                                   )
