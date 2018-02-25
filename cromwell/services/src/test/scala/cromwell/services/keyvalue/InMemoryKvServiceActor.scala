package cromwell.services.keyvalue

import cromwell.services.keyvalue.KeyValueServiceActor._

import scala.concurrent.{ExecutionContextExecutor, Future}

final class InMemoryKvServiceActor extends KeyValueServiceActor {
  override implicit val ec: ExecutionContextExecutor = context.dispatcher

  var kvStore = Map.empty[ScopedKey, Option[String]]

  override def doGet(get: KvGet): Future[KvResponse] = kvStore.get(get.key).map(KvPair(get.key, _)) match {
    case Some(kvPair) => Future.successful(kvPair)
    case None => Future.successful(KvKeyLookupFailed(get))
  }

  override def doPut(put: KvPut): Future[KvResponse] = {
    kvStore += (put.key -> put.pair.value)
    Future.successful(KvPutSuccess(put))
  }
}