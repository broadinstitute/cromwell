package cromwell.engine.workflow.tokens

import cromwell.backend.standard.GroupMetricsActor
import cromwell.backend.standard.GroupMetricsActor.{GetQuotaExhaustedGroups, GetQuotaExhaustedGroupsSuccess}
import cromwell.services.EngineServicesStore.engineDatabaseInterface

object TokenDispenserUtils {

  class TestGroupMetricsActor extends GroupMetricsActor(engineDatabaseInterface, 15) {
    override def receive: Receive = { case GetQuotaExhaustedGroups =>
      sender() ! GetQuotaExhaustedGroupsSuccess(List.empty)
    }
  }
}
