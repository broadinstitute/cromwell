package cromwell.engine.db

import cromwell.binding.FullyQualifiedName

trait CallBackendInfo {
  val callFqn: FullyQualifiedName
  val status: CallStatus
}
