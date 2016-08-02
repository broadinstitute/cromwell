package cromwell.services.keyvalue.impl

import scala.concurrent.{Future, ExecutionContext}
import akka.actor.ActorSystem

import cromwell.database.Database
import cromwell.util.DatabaseUtil._

import cromwell.core.{JobKey, WorkflowId}
import cromwell.core.ExecutionIndex._
import wdl4s.Scope

trait BackendKeyValueDatabaseAccess { this: Database =>

  def getBackendValueByKey(workflowId: WorkflowId, call: Scope, callIndex: Option[Int], attempt: Int, key: String)
                           (implicit ec: ExecutionContext): Future[Option[String]] = {
    databaseInterface.queryBackendStoreValueByStoreKey(workflowId.toString, call.fullyQualifiedName, callIndex.fromIndex, attempt, key)
  }

  def updateBackendKeyValuePair(workflowId: WorkflowId,
                                callKey: JobKey,
                                backendStoreKey: String,
                                backendStoreValue: String)(implicit ec: ExecutionContext, actorSystem: ActorSystem): Future[Unit] = {

    withRetry (() =>
      databaseInterface.addBackendStoreKeyValuePair(workflowId.toString, callKey.scope.fullyQualifiedName, callKey.index.fromIndex, callKey.attempt, backendStoreKey, backendStoreValue)
    )
  }

}
