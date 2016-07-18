package cromwell

import java.time.OffsetDateTime

import cromwell.core.JobOutput
import wdl4s._
import wdl4s.values.WdlValue

import scala.language.implicitConversions
import scala.util.{Failure, Success, Try}

package object engine {

  final case class AbortFunction(function: () => Unit)
  final case class AbortRegistrationFunction(register: AbortFunction => Unit)


  final case class FailureEventEntry(failure: String, timestamp: OffsetDateTime)
  final case class CallAttempt(fqn: FullyQualifiedName, attempt: Int)

  implicit class EnhancedFullyQualifiedName(val fqn: FullyQualifiedName) extends AnyVal {
    def scopeAndVariableName: (String, String) = {
      val array = fqn.split("\\.(?=[^\\.]+$)")
      (array(0), array(1))
    }
  }

  implicit class EnhancedCallOutputMap[A](val m: Map[A, JobOutput]) extends AnyVal {
    def mapToValues: Map[A, WdlValue] = m map {
      case (k, JobOutput(wdlValue)) => (k, wdlValue)
    }
  }

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



}
