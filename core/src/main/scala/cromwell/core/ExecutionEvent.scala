package cromwell.core

import java.time.OffsetDateTime

final case class ExecutionEvent(name: String, offsetDateTime: OffsetDateTime)

object ExecutionEvent {
  def apply(name: String): ExecutionEvent = ExecutionEvent(name, OffsetDateTime.now())
}
