package cromwell.backend.google.batch.models

import cromwell.backend.BackendJobDescriptor
import cromwell.core.path.Path

case class CreateGcpBatchParameters(jobDescriptor: BackendJobDescriptor,
                                    runtimeAttributes: GcpBatchRuntimeAttributes,
                                    batchAttributes: GcpBatchConfigurationAttributes,
                                    projectId: String,
                                    region: String,
                                    logfile: Path
)
