package cromwell.backend.google.batch.models

import cromwell.backend.google.batch.api.GcpBatchRequestFactory.CreateBatchJobParameters
import cromwell.core.WorkflowId

case class GcpBatchRequest(workflowId: WorkflowId,
                           createParameters: CreateBatchJobParameters,
                           jobName: String,
                           gcpBatchParameters: CreateGcpBatchParameters
)
