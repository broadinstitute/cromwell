package cromwell.services.keyvalue

import akka.actor.{Actor, ActorRef, Props}
import com.typesafe.config.Config
import cromwell.core.{WorkflowId, JobKey}
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage

import scala.concurrent.Future
import scala.util.{Failure, Success}

object KeyValueService {
  trait KvMessage
  trait KvAction extends KvMessage with ServiceRegistryMessage {
    def serviceName = "KeyValue"
  }
  trait KvResponse extends KvMessage

  case class ScopedKey(workflowId: WorkflowId, jobKey: JobKey, key: String)
  case class KvPut(pair: KvPair) extends KvAction
  case class KvGet(key: ScopedKey) extends KvAction

  case class KvPair(key: ScopedKey, value: Option[String]) extends KvResponse
  case class KvFailure(action: KvAction, failure: Throwable) extends KvResponse
  case class KvKeyLookupFailed(action: KvGet) extends KvResponse
  case class KvPutSuccess(action: KvPut) extends KvResponse
}