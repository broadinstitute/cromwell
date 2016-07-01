package cromwell.backend.impl.htcondor.caching.model

import cromwell.backend.BackendJobExecutionActor.SucceededResponse

case class CachedExecutionResult(hash: String, succeededResponse: SucceededResponse)

