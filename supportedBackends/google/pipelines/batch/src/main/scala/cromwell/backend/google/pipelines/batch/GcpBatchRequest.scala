package cromwell.backend.google.pipelines.batch

import cromwell.backend.google.pipelines.batch.api.GcpBatchRequestFactory.CreatePipelineParameters
import cromwell.core.WorkflowId

case class GcpBatchRequest(workflowId: WorkflowId,
                           createParameters: CreatePipelineParameters,
                           jobName: String,
                           gcpBatchCommand: String,
                           gcpBatchParameters: CreateGcpBatchParameters)
