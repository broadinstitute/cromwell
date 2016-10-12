package cromwell.backend.impl.htcondor.caching.model

import cromwell.backend.BackendJobExecutionActor.JobSucceededResponse

case class CachedExecutionResult(hash: String, succeededResponse: JobSucceededResponse)

