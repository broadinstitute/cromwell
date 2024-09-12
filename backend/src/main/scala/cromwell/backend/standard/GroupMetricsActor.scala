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
import scala.concurrent.Future
import scala.concurrent.duration.FiniteDuration
import scala.util.{Failure, Success}

class GroupMetricsActor(engineDbInterface: EngineSqlDatabase, quotaExhaustionThresholdInMins: Long, loggingInterval: FiniteDuration)
    extends Actor
    with ActorLogging {

  implicit val ec: MessageDispatcher = context.system.dispatchers.lookup(Dispatcher.EngineDispatcher)

  // initial schedule for logging exhausted groups
  context.system.scheduler.scheduleOnce(loggingInterval)(self ! LogQuotaExhaustedGroups)
  log.info(s"${this.getClass.getSimpleName} configured to log groups experiencing quota exhaustion at interval ${loggingInterval.toString()}.")

  override def receive: Receive = {
    case RecordGroupQuotaExhaustion(group) =>
      val groupMetricsEntry = GroupMetricsEntry(group, OffsetDateTime.now.toSystemTimestamp)
      engineDbInterface.recordGroupMetricsEntry(groupMetricsEntry)
      ()
    case GetQuotaExhaustedGroups =>
      val respondTo: ActorRef = sender()
      getQuotaExhaustedGroups() onComplete {
        case Success(quotaExhaustedGroups) => respondTo ! GetQuotaExhaustedGroupsSuccess(quotaExhaustedGroups.toList)
        case Failure(exception) => respondTo ! GetQuotaExhaustedGroupsFailure(exception.getMessage)
      }
    case LogQuotaExhaustedGroups => getQuotaExhaustedGroups() onComplete {
        case Success(quotaExhaustedGroups) =>
          log.info(s"Hog groups currently experiencing quota exhaustion: ${quotaExhaustedGroups.length}. Group IDs: ${quotaExhaustedGroups.toList.toString()}")
        case Failure(exception) =>
          log.info(s"Something went wrong when fetching quota exhausted groups for logging. Will retry in ${loggingInterval.toString()}. Exception: ${exception.getMessage}")
      }
      // schedule next logging
      context.system.scheduler.scheduleOnce(loggingInterval)(self ! LogQuotaExhaustedGroups)
      ()
    case other =>
      log.error(
        s"Programmer Error: Unexpected message ${other.toPrettyElidedString(1000)} received by ${this.self.path.name}."
      )
  }

  private def getQuotaExhaustedGroups(): Future[Seq[String]] = {
    // for a group in the GROUP_METRICS_ENTRY table, if the 'quota_exhaustion_detected' timestamp hasn't
    // been updated in last X minutes it is no longer experiencing cloud quota exhaustion
    val currentTimestampMinusDelay = OffsetDateTime.now().minusMinutes(quotaExhaustionThresholdInMins)
    engineDbInterface.getQuotaExhaustedGroups(currentTimestampMinusDelay.toSystemTimestamp)
  }
}

object GroupMetricsActor {

  // Requests
  sealed trait GroupMetricsActorMessage
  case class RecordGroupQuotaExhaustion(group: String) extends GroupMetricsActorMessage
  case object GetQuotaExhaustedGroups extends GroupMetricsActorMessage
  case object LogQuotaExhaustedGroups extends GroupMetricsActorMessage

  // Responses
  sealed trait GetQuotaExhaustedGroupsResponse
  case class GetQuotaExhaustedGroupsSuccess(quotaExhaustedGroups: List[String]) extends GetQuotaExhaustedGroupsResponse
  case class GetQuotaExhaustedGroupsFailure(errorMsg: String) extends GetQuotaExhaustedGroupsResponse

  def props(engineDbInterface: EngineSqlDatabase, quotaExhaustionThresholdInMins: Long, loggingInterval: FiniteDuration): Props =
    Props(new GroupMetricsActor(engineDbInterface, quotaExhaustionThresholdInMins, loggingInterval)).withDispatcher(EngineDispatcher)
}
