package cromwell.services

import akka.actor.{Actor, Props}
import akka.pattern.pipe
import com.typesafe.config.Config
import cromwell.engine.WorkflowId
import cromwell.engine.db.DataAccess
import cromwell.engine.workflow.BackendCallKey
import cromwell.services.KeyValueServiceActor._

import scala.concurrent.Future

object KeyValueServiceActor {
  trait KeyValueMessage
  trait KeyValueAction extends KeyValueMessage
  trait KeyValueResponse extends KeyValueMessage

  case class Put(backendCallKey: BackendCallKey, key: String, value: String) extends KeyValueAction
  case class Get(backendCallKey: BackendCallKey, key: String) extends KeyValueAction

  case class KeyValuePair(key: String, value: Option[String]) extends KeyValueResponse
  case class KeyValueLookupFailed(failedAction: Get, failure: Throwable) extends KeyValueResponse

  def props(dataAccess: DataAccess, workflowId: WorkflowId, serviceConfigPath: String, config: Config) = {
    Props(KeyValueServiceActor(dataAccess, workflowId, serviceConfigPath, config))
  }
}

case class KeyValueServiceActor(dataAccess: DataAccess, workflowId: WorkflowId, serviceConfigPath: String, entireConfig: Config) extends Actor {
  private implicit val ec = context.dispatcher

  def receive = {
    case action: Put =>
      println("put action")
      pipe(put(action)) to sender
    case action: Get => pipe(get(action)) to sender
  }

  private def put(put: Put): Future[KeyValueResponse] = {
    dataAccess.updateExecutionInfo(workflowId, put.backendCallKey, put.key, Option(put.value)).map(
      _ => {
        println("put complete")
        KeyValuePair(put.key, Option(put.value))
      }
    )
  }

  private def get(get: Get): Future[KeyValueResponse] = {
    dataAccess.getExecutionInfos(workflowId, get.backendCallKey.scope, get.backendCallKey.attempt) map {
      case infos => infos.find(_.key == get.key) match {
        case Some(info) => KeyValuePair(get.key, info.value)
        case None => KeyValueLookupFailed(get, new Throwable(s"Cannot lookup key '${get.key}' from $workflowId:${get.backendCallKey.tag}"))
      }
    }
  }
}
