package cromwell.services.metadata.impl

import akka.actor.{ActorLogging, ActorRef, Props}
import cats.data.NonEmptyVector
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.core.Mailbox.PriorityMailbox
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.services.metadata.MetadataEvent
import cromwell.services.metadata.MetadataService._
import cromwell.services.{EnhancedBatchActor, MetadataServicesStore}

import scala.concurrent.duration._
import scala.util.{Failure, Success}


class WriteMetadataActor(override val batchSize: Int,
                         override val flushRate: FiniteDuration,
                         override val serviceRegistryActor: ActorRef,
                         override val threshold: Int)
  extends EnhancedBatchActor[MetadataWriteAction](flushRate, batchSize)
    with ActorLogging
    with MetadataDatabaseAccess
    with MetadataServicesStore {

  override def process(e: NonEmptyVector[MetadataWriteAction]) = instrumentedProcess {
    val empty = (Vector.empty[MetadataEvent], List.empty[(Iterable[MetadataEvent], ActorRef)])

    val (putWithoutResponse, putWithResponse) = e.foldLeft(empty)({
      case ((putEvents, putAndRespondEvents), action: PutMetadataAction) =>
        (putEvents ++ action.events, putAndRespondEvents)
      case ((putEvents, putAndRespondEvents), action: PutMetadataActionAndRespond) =>
        (putEvents, putAndRespondEvents :+ (action.events -> action.replyTo))
    })
    val allPutEvents: Iterable[MetadataEvent] = putWithoutResponse ++ putWithResponse.flatMap(_._1)
    val dbAction = addMetadataEvents(allPutEvents)

    dbAction onComplete {
      case Success(_) =>
        putWithResponse foreach { case (ev, replyTo) => replyTo ! MetadataWriteSuccess(ev) }
      case Failure(regerts) =>

        val (outOfLives, stillGood) = e.toVector.partition(_.ttl <= 1)

        handleOutOfLives(outOfLives, regerts)
        handleEventsToReconsider(stillGood)
    }

    dbAction.map(_ => allPutEvents.size)
  }

  private def handleOutOfLives(writeActions: Vector[MetadataWriteAction], reason: Throwable): Unit = {
    val workflowMetadataFailureCounts = writeActions.toVector.flatMap(_.events).groupBy(x => x.key.workflowId).map { case (wfid, list) => s"$wfid: ${list.size}" }
    log.error("Metadata events have been permanently dropped for the following workflows: " + workflowMetadataFailureCounts.mkString(","))

    writeActions foreach {
      case PutMetadataActionAndRespond(ev, replyTo, _) => replyTo ! MetadataWriteFailure(reason, ev)
      case _: PutMetadataAction => () // We need to satisfy the exhaustive match but there's nothing special to do here
    }
  }

  private def handleEventsToReconsider(writeActions: Vector[MetadataWriteAction]): Unit = {
    val workflowMetadataFailureCounts = writeActions.toVector.flatMap(_.events).groupBy(x => x.key.workflowId).map { case (wfid, list) => s"$wfid: ${list.size}" }
    log.error("Metadata event writes failed for the following workflows but will be reconsidered: " + workflowMetadataFailureCounts.mkString(","))

    writeActions foreach {
      case action: PutMetadataAction => self ! action.copy(ttl = action.ttl - 1)
      case action: PutMetadataActionAndRespond => self ! action.copy(ttl = action.ttl - 1)
    }
  }

  // EnhancedBatchActor overrides
  override def receive = enhancedReceive.orElse(super.receive)
  override protected def weightFunction(command: MetadataWriteAction) = command.size
  override protected def instrumentationPath = MetadataServiceActor.MetadataInstrumentationPrefix
  override protected def instrumentationPrefix = InstrumentationPrefixes.ServicesPrefix
  def commandToData(snd: ActorRef): PartialFunction[Any, MetadataWriteAction] = {
    case command: MetadataWriteAction => command
  }
}

object WriteMetadataActor {
  val MetadataWritePath = MetadataServiceActor.MetadataInstrumentationPrefix.::("writes")

  def props(dbBatchSize: Int,
            flushRate: FiniteDuration,
            serviceRegistryActor: ActorRef,
            threshold: Int): Props =
    Props(new WriteMetadataActor(dbBatchSize, flushRate, serviceRegistryActor, threshold))
      .withDispatcher(ServiceDispatcher)
      .withMailbox(PriorityMailbox)
}
