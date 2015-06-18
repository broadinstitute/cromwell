package cromwell.engine.db

import cromwell.binding.FullyQualifiedName

case class LocalCallInfo(callFqn: FullyQualifiedName, status: CallStatus,
                         processId: Int, command: String, resultCode: Int) extends CallInfo
