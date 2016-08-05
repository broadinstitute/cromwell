package cromwell.services.keyvalue.impl

import akka.actor.Props
import com.typesafe.config.Config
import cromwell.database.CromwellDatabase
import cromwell.services.keyvalue.KeyValueServiceActor
import cromwell.services.keyvalue.KeyValueServiceActor._

import scala.concurrent.Future

object SqlKeyValueServiceActor {
  def props(serviceConfig: Config, globalConfig: Config) = Props(SqlKeyValueServiceActor(serviceConfig, globalConfig))
}

case class SqlKeyValueServiceActor(override val serviceConfig: Config, override val globalConfig: Config) extends KeyValueServiceActor with BackendKeyValueDatabaseAccess with CromwellDatabase {
  override implicit val ec = context.dispatcher
  private implicit val system = context.system

  override def doPut(put: KvPut): Future[KvResponse] = {
    put.pair.value match {
      case Some(backendVal) => updateBackendKeyValuePair(put.pair.key.workflowId,
        put.pair.key.jobKey,
        put.pair.key.key,
        put.pair.value.get).map(_ => KvPutSuccess(put))
      case None => Future.successful(KvFailure(put, new RuntimeException(s"Failed to find the value associated to key: ${put.pair.key.key}. This key cannot be added to the BackendKVStore.")))
    }
  }

  override def doGet(get: KvGet): Future[KvResponse] = {
    val backendValue = getBackendValueByKey(
      get.key.workflowId,
      get.key.jobKey.scope,
      get.key.jobKey.index,
      get.key.jobKey.attempt,
      get.key.key
    )

    backendValue map {
      case Some(maybeValue) => KvPair(get.key, Option(maybeValue))
      case None => KvKeyLookupFailed(get)
    }
  }
}
