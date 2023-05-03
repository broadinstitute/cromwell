package cromwell.backend.google.batch.models

import cromwell.backend.google.batch.api.GcpBatchRequestFactory.CreatePipelineParameters
import cromwell.core.WorkflowId

case class GcpBatchRequest(workflowId: WorkflowId,
                           createParameters: CreatePipelineParameters,
                           jobName: String,
                           gcpBatchParameters: CreateGcpBatchParameters)
