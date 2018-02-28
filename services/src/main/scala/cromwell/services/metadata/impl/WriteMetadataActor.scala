package cromwell.services.metadata.impl

import akka.actor.{ActorLogging, ActorRef, Props}
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.core.Mailbox.PriorityMailbox
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.services.MetadataServicesStore
import cromwell.services.instrumentation.{CromwellInstrumentationActor, InstrumentedBatchActor}
import cromwell.services.metadata.MetadataEvent
import cromwell.services.metadata.MetadataService._

import scala.concurrent.duration._
import scala.util.{Failure, Success}


class WriteMetadataActor(override val batchSize: Int,
                         override val flushRate: FiniteDuration,
                         override val serviceRegistryActor: ActorRef)
  extends InstrumentedBatchActor[MetadataWriteAction](flushRate, batchSize,
    MetadataServiceActor.MetadataInstrumentationPrefix, InstrumentationPrefixes.ServicesPrefix) with ActorLogging with
    MetadataDatabaseAccess with MetadataServicesStore with CromwellInstrumentationActor {

  def commandToData(snd: ActorRef): PartialFunction[Any, MetadataWriteAction] = {
    case command: MetadataWriteAction => command
  }

  override def processInner(e: Vector[MetadataWriteAction]) = {
    val empty = (Vector.empty[MetadataEvent], Map.empty[Iterable[MetadataEvent], ActorRef])

    val (putWithoutResponse, putWithResponse) = e.foldLeft(empty)({
      case ((putEvents, putAndRespondEvents), action: PutMetadataAction) =>
        (putEvents ++ action.events, putAndRespondEvents)
      case ((putEvents, putAndRespondEvents), action: PutMetadataActionAndRespond) =>
        (putEvents, putAndRespondEvents + (action.events -> action.replyTo))
    })
    val allPutEvents: Iterable[MetadataEvent] = putWithoutResponse ++ putWithResponse.keys.flatten
    val dbAction = addMetadataEvents(allPutEvents)

    dbAction onComplete {
      case Success(_) =>
        putWithResponse foreach { case (ev, replyTo) => replyTo ! MetadataWriteSuccess(ev) }
      case Failure(regerts) =>
        putWithResponse foreach { case (ev, replyTo) => replyTo ! MetadataWriteFailure(regerts, ev) }
    }

    dbAction.map(_ => allPutEvents.size)
  }

  override protected def weightFunction(command: MetadataWriteAction) = command.size
}

object WriteMetadataActor {
  val MetadataWritePath = MetadataServiceActor.MetadataInstrumentationPrefix.::("writes")

  def props(dbBatchSize: Int,
            flushRate: FiniteDuration,
            serviceRegistryActor: ActorRef): Props =
    Props(new WriteMetadataActor(dbBatchSize, flushRate, serviceRegistryActor))
      .withDispatcher(ServiceDispatcher)
      .withMailbox(PriorityMailbox)
}
