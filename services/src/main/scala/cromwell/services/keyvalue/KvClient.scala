package cromwell.services.keyvalue

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.services.keyvalue.KeyValueServiceActor._

import scala.concurrent.{ExecutionContext, Future, Promise}

trait KvClient { this: Actor with ActorLogging =>

  def serviceRegistryActor: ActorRef
  private[keyvalue] var currentKvClientRequests: Map[ScopedKey, Promise[KvResponse]] = Map.empty

  final def makeKvRequest(actions: Seq[KvAction])(implicit ec: ExecutionContext): Future[Seq[KvResponse]] = {
    if (actions.exists(action => currentKvClientRequests.contains(action.key))) {
      val msg = "Programmer Error! KvClient does not support multiple KvActions active for the same ScopedKey concurrently. Mi Scusi!"
      log.error(msg)
      Future.failed(new RuntimeException(msg))
    } else {
      createResponseSet(actions)
    }
  }

  final def kvClientReceive: Actor.Receive = {
    case response: KvResponse => fulfillOrLog(response)
  }

  private def createResponseSet(newActions: Seq[KvAction])(implicit ec: ExecutionContext) = {
    val actionsAndPromises = newActions.map(a => a.key -> Promise[KvResponse]())
    currentKvClientRequests ++= actionsAndPromises.toMap
    newActions foreach { serviceRegistryActor ! _ }
    Future.sequence(actionsAndPromises.map(_._2.future))
  }

  private def fulfillOrLog(response: KvResponse) = currentKvClientRequests.get(response.key) match {
    case Some(fulfilledPromise) =>
      fulfilledPromise.success(response)
      currentKvClientRequests -= response.key
    case None => log.error(s"Programmer Error: Got a KV response for a request that was never sent: $response. Did you use the KV store without KvClient? Current key set: ${currentKvClientRequests.keys.mkString("")}")
  }
}

