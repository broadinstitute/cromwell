package cromwell.backend.google.pipelines.batch.models

import cromwell.backend.google.pipelines.batch.api.GcpBatchRequestFactory.CreatePipelineParameters
import cromwell.core.WorkflowId

case class GcpBatchRequest(workflowId: WorkflowId,
                           createParameters: CreatePipelineParameters,
                           jobName: String,
                           gcpBatchParameters: CreateGcpBatchParameters)
