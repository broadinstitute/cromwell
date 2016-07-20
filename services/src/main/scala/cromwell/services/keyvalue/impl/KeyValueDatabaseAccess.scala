package cromwell.services.keyvalue.impl

import akka.actor.ActorSystem
import cromwell.core.ExecutionIndex._
import cromwell.util.DatabaseUtil._
import cromwell.core.{JobKey, WorkflowId}
import cromwell.database.Database
import wdl4s.Scope

import scala.concurrent.{ExecutionContext, Future}

trait KeyValueDatabaseAccess { this: Database =>

  def getExecutionInfoByKey(workflowId: WorkflowId, call: Scope, attempt: Int, key: String)
                           (implicit ec: ExecutionContext): Future[Option[Option[String]]] = {
    databaseInterface.getExecutionInfoByKey(workflowId.toString, call.fullyQualifiedName, attempt, key)
  }

  /**
    * TODO: the interface for retrying SQL commands that might fail because
    * of transient reasons (e.g. upsertExecutionInfo and upsertRuntimeAttributes)
    * could be made better.  Also it'd be nice if it were transparently
    * turned on/off for all methods in dataAccess.
    *
    * https://github.com/broadinstitute/cromwell/issues/693
    */
  def upsertExecutionInfo(workflowId: WorkflowId,
                          callKey: JobKey,
                          keyValues: Map[String, Option[String]])(implicit ec: ExecutionContext, actorSystem: ActorSystem): Future[Unit] = {
    withRetry {
      databaseInterface.upsertExecutionInfo(workflowId.toString, callKey.scope.fullyQualifiedName, callKey.index.fromIndex, callKey.attempt, keyValues)
    }
  }
}
