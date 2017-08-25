package cromwell.services.keyvalue.impl

import akka.actor.ActorSystem
import cromwell.core.ExecutionIndex._
import cromwell.core.WorkflowId
import cromwell.database.sql.tables.JobKeyValueEntry
import cromwell.services.EngineServicesStore
import cromwell.services.keyvalue.KeyValueServiceActor.KvJobKey
import cromwell.util.DatabaseUtil._

import scala.concurrent.{ExecutionContext, Future}

trait BackendKeyValueDatabaseAccess {

  def getBackendValueByKey(workflowId: WorkflowId, jobKey: KvJobKey, key: String)
                           (implicit ec: ExecutionContext): Future[Option[String]] = {
    EngineServicesStore.engineDatabaseInterface.queryStoreValue(
      workflowId.toString, jobKey.callFqn, jobKey.callIndex.fromIndex, jobKey.callAttempt, key)
  }

  def updateBackendKeyValuePair(workflowId: WorkflowId,
                                jobKey: KvJobKey,
                                backendStoreKey: String,
                                backendStoreValue: String)(implicit ec: ExecutionContext, actorSystem: ActorSystem): Future[Unit] = {
    val jobKeyValueEntry = JobKeyValueEntry(workflowId.toString, jobKey.callFqn, jobKey.callIndex.fromIndex,
      jobKey.callAttempt, backendStoreKey, backendStoreValue)
    withRetry(() => EngineServicesStore.engineDatabaseInterface.addJobKeyValueEntry(jobKeyValueEntry))
  }
}
