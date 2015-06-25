package cromwell.engine.db

import cromwell.binding.FullyQualifiedName

case class LocalCallBackendInfo(callFqn: FullyQualifiedName, status: CallStatus,
                         processId: Int, command: String, resultCode: Int) extends CallBackendInfo
