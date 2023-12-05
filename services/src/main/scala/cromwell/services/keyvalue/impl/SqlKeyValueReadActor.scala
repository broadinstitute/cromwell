package cromwell.services.keyvalue.impl

import akka.actor.{ActorRef, Props}
import cromwell.services.keyvalue.KeyValueServiceActor.{KvKeyLookupFailed, KvPair}
import cromwell.services.keyvalue.{KeyValueReadActor, KeyValueServiceActor}

object SqlKeyValueReadActor {
  def props(threshold: Int, serviceRegistryActor: ActorRef) =
    Props(new SqlKeyValueReadActor(threshold, serviceRegistryActor))
}

class SqlKeyValueReadActor(threshold: Int, serviceRegistryActor: ActorRef)
    extends KeyValueReadActor(threshold, serviceRegistryActor)
    with BackendKeyValueDatabaseAccess {
  override def processGet(get: KeyValueServiceActor.KvGet) = {
    val backendValue = getBackendValueByKey(
      get.key.workflowId,
      get.key.jobKey,
      get.key.key
    )

    backendValue map {
      case Some(maybeValue) => KvPair(get.key, maybeValue)
      case None => KvKeyLookupFailed(get)
    }
  }
}
