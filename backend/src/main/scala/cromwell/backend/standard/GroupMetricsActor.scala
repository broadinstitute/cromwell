package cromwell.backend.standard

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.dispatch.MessageDispatcher
import common.util.StringUtil.EnhancedToStringable
import cromwell.backend.standard.GroupMetricsActor._
import cromwell.core.Dispatcher
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.database.sql.EngineSqlDatabase
import cromwell.database.sql.SqlConverters.OffsetDateTimeToSystemTimestamp
import cromwell.database.sql.tables.GroupMetricsEntry

import java.time.OffsetDateTime
import scala.util.{Failure, Success}

class GroupMetricsActor(engineDbInterface: EngineSqlDatabase) extends Actor with ActorLogging {

  implicit val ec: MessageDispatcher = context.system.dispatchers.lookup(Dispatcher.EngineDispatcher)

  final private val QUOTA_EXHAUSTION_THRESHOLD_IN_SECS = 15 * 60 // 15 minutes

  override def receive: Receive = {
    case RecordGroupQuotaExhaustion(group) =>
      val groupMetricsEntry = GroupMetricsEntry(group, OffsetDateTime.now.toSystemTimestamp)
      engineDbInterface.recordGroupMetricsEntry(groupMetricsEntry)
      ()
    case GetQuotaExhaustedGroups =>
      val respondTo: ActorRef = sender()

      // for a group in the GROUP_METRICS_ENTRY table, if the 'quota_exhaustion_detected' timestamp hasn't
      // been updated in last 15 minutes it is no longer experiencing cloud quota exhaustion
      val currentTimestampMinusDelay = OffsetDateTime.now().minusSeconds(QUOTA_EXHAUSTION_THRESHOLD_IN_SECS)
      engineDbInterface.getQuotaExhaustedGroups(currentTimestampMinusDelay.toSystemTimestamp) onComplete {
        case Success(quotaExhaustedGroups) => respondTo ! GetQuotaExhaustedGroupsSuccess(quotaExhaustedGroups.toList)
        case Failure(exception) => respondTo ! GetQuotaExhaustedGroupsFailure(exception.getMessage)
      }
    case other =>
      log.error(
        s"Programmer Error: Unexpected message ${other.toPrettyElidedString(1000)} received by ${this.self.path.name}."
      )
  }
}

object GroupMetricsActor {

  // Requests
  sealed trait GroupMetricsActorMessage
  case class RecordGroupQuotaExhaustion(group: String) extends GroupMetricsActorMessage
  case object GetQuotaExhaustedGroups extends GroupMetricsActorMessage

  // Responses
  sealed trait GetQuotaExhaustedGroupsResponse
  case class GetQuotaExhaustedGroupsSuccess(quotaExhaustedGroups: List[String]) extends GetQuotaExhaustedGroupsResponse
  case class GetQuotaExhaustedGroupsFailure(errorMsg: String) extends GetQuotaExhaustedGroupsResponse

  def props(engineDbInterface: EngineSqlDatabase): Props =
    Props(new GroupMetricsActor(engineDbInterface)).withDispatcher(EngineDispatcher)
}
