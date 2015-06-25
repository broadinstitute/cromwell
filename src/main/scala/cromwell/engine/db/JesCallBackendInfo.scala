package cromwell.engine.db

import cromwell.binding.FullyQualifiedName

// Assumes that `jes_job` or `local_job` are 1-1 with `execution`
case class JesCallBackendInfo(callFqn: FullyQualifiedName, status: CallStatus,
                       jesId: JesId, jesStatus: JesStatus) extends CallBackendInfo
