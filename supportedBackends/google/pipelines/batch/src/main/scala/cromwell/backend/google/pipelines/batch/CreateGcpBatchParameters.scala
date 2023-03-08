package cromwell.backend.google.pipelines.batch

import cromwell.backend.BackendJobDescriptor

case class CreateGcpBatchParameters(jobDescriptor: BackendJobDescriptor,
                                    runtimeAttributes: GcpBatchRuntimeAttributes,
                                    dockerImage: String,
                                    projectId: String,
                                    region: String
                                   )
