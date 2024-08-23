package cromwell.backend.standard

import akka.actor.{Actor, ActorLogging, Props}
import akka.dispatch.MessageDispatcher
import common.util.StringUtil.EnhancedToStringable
import cromwell.backend.standard.GroupMetricsActor.RecordGroupQuotaExhaustion
import cromwell.core.Dispatcher
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.database.sql.EngineSqlDatabase
import cromwell.database.sql.SqlConverters.OffsetDateTimeToSystemTimestamp
import cromwell.database.sql.tables.GroupMetricsEntry

import java.time.OffsetDateTime

class GroupMetricsActor(engineDbInterface: EngineSqlDatabase) extends Actor with ActorLogging {

  implicit val ec: MessageDispatcher = context.system.dispatchers.lookup(Dispatcher.EngineDispatcher)

  override def receive: Receive = {
    case RecordGroupQuotaExhaustion(group) =>
      val groupMetricsEntry = GroupMetricsEntry(group, OffsetDateTime.now.toSystemTimestamp)
      engineDbInterface.recordGroupMetricsEntry(groupMetricsEntry)
      ()
    case other =>
      log.error(
        s"Programmer Error: Unexpected message ${other.toPrettyElidedString(1000)} received by ${this.self.path.name}."
      )
  }
}

object GroupMetricsActor {

  sealed trait GroupMetricsActorMessage

  case class RecordGroupQuotaExhaustion(group: String) extends GroupMetricsActorMessage

  def props(engineDbInterface: EngineSqlDatabase): Props =
    Props(new GroupMetricsActor(engineDbInterface)).withDispatcher(EngineDispatcher)
}
