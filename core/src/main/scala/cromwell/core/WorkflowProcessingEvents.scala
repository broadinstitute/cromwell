package cromwell.core

import java.time.OffsetDateTime

object WorkflowProcessingEvents {

  sealed trait WithSimpleCaseObjectName {
    def simpleCaseObjectName: String = {
      val x = getClass.getName.dropRight(1)
      x.substring(x.lastIndexOf('$') + 1)
    }
  }

  object EventKey {
    sealed trait Key extends WithSimpleCaseObjectName {
      def key: String = {
        val upper = simpleCaseObjectName
        upper(0).toLower +: upper.tail
      }
    }
    case object Description extends Key
    case object CromwellId extends Key
    case object Timestamp extends Key
    case object CromwellVersion extends Key
  }

  object DescriptionEventValue {
    sealed trait Value extends WithSimpleCaseObjectName {
      def value: String = simpleCaseObjectName
    }
    case object PickedUp extends Value
    case object Released extends Value
    case object Finished extends Value
  }

  val ProcessingEventsKey = "workflowProcessingEvents"
}

case class WorkflowProcessingEvent(cromwellId: String, description: String, timestamp: OffsetDateTime, cromwellVersion: String)
