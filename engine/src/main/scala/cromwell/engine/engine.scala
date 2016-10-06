package cromwell.engine

import java.time.OffsetDateTime

import wdl4s._

import scala.util.{Failure, Success, Try}

final case class AbortFunction(function: () => Unit)
final case class AbortRegistrationFunction(register: AbortFunction => Unit)

final case class FailureEventEntry(failure: String, timestamp: OffsetDateTime)
final case class CallAttempt(fqn: FullyQualifiedName, attempt: Int)

object WorkflowFailureMode {
  def tryParse(mode: String): Try[WorkflowFailureMode] = {
    val modes = Seq(ContinueWhilePossible, NoNewCalls)
    modes find { _.toString.equalsIgnoreCase(mode) } map { Success(_) } getOrElse Failure(new Exception(s"Invalid workflow failure mode: $mode"))
  }
}
sealed trait WorkflowFailureMode {
  def allowNewCallsAfterFailure: Boolean
}
case object ContinueWhilePossible extends WorkflowFailureMode { override val allowNewCallsAfterFailure = true }
case object NoNewCalls extends WorkflowFailureMode { override val allowNewCallsAfterFailure = false }