package cromwell.engine.db

import cromwell.binding.FullyQualifiedName

trait CallInfo {
  val callFqn: FullyQualifiedName
  val status: CallStatus
}
