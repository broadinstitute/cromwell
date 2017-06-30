package cromwell.services.metadata.impl

import akka.actor.{ActorLogging, ActorRef, LoggingFSM, Props}
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.core.actor.BatchingDbWriter
import cromwell.core.actor.BatchingDbWriter._
import cromwell.services.SingletonServicesStore
import cromwell.services.metadata.MetadataEvent
import cromwell.services.metadata.MetadataService._

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.util.{Failure, Success}


class WriteMetadataActor(batchSize: Int, flushRate: FiniteDuration)
  extends LoggingFSM[BatchingDbWriterState, BatchingDbWriter.BatchingDbWriterData] with ActorLogging with
  MetadataDatabaseAccess with SingletonServicesStore {
  import WriteMetadataActor._

  implicit val ec: ExecutionContext = context.dispatcher

  override def preStart(): Unit = {
    context.system.scheduler.schedule(0.seconds, flushRate, self, ScheduledFlushToDb)
    super.preStart()
  }


  startWith(WaitingToWrite, NoData)

  when(WaitingToWrite) {
    case Event(PutMetadataAction(events), curData) =>
      curData.addData(events) match {
        case newData: HasData[_] if newData.length > batchSize => goto(WritingToDb) using newData
        case newData => stay using newData
      }
    case Event(ScheduledFlushToDb, curData) =>
      log.debug("Initiating periodic metadata flush to DB")
      goto(WritingToDb) using curData
    case Event(CheckPendingWrites, NoData) =>
      sender() ! NoPendingWrites
      stay()
    case Event(CheckPendingWrites, _: HasData[_]) =>
      sender() ! HasPendingWrites
      stay()
    case Event(e @ PutMetadataActionAndRespond(events, replyTo), curData) =>
      curData.addData(e) match {
        case newData: HasData[_] if newData.length > batchSize => goto(WritingToDb) using newData
        case newData => stay using newData
      }
  }

  when(WritingToDb) {
    case Event(CheckPendingWrites, _) => 
      sender() ! HasPendingWrites
      stay()
    case Event(ScheduledFlushToDb, curData) => stay using curData
    case Event(PutMetadataAction(events), curData) => stay using curData.addData(events)
    case Event(FlushBatchToDb, NoData) =>
      log.debug("Attempted metadata flush to DB but had nothing to write")
      goto(WaitingToWrite) using NoData
    case Event(FlushBatchToDb, HasData(e)) =>
      log.debug("Flushing {} metadata events to the DB", e.length)
      // blech
      //Partitioning the current data into put events that require a response and those that don't
      val empty = (List.empty[MetadataEvent], Map.empty[Iterable[MetadataEvent], ActorRef])
      val (putWithoutResponse, putWithResponse) = e.toVector.foldLeft(empty)({
        case ((putEvents, putAndRespondEvents), events) =>
          events match {
            case putEvent: MetadataEvent => (putEvents :+ putEvent, putAndRespondEvents)
            case PutMetadataActionAndRespond(ev, replyTo) => (putEvents, putAndRespondEvents + (ev -> replyTo))
          }
      })
      val allPutEvents: Iterable[MetadataEvent] = putWithoutResponse ++ putWithResponse.keys.flatten
      addMetadataEvents(allPutEvents) onComplete {
        case Success(_) =>
          self ! DbWriteComplete
          putWithResponse foreach { case(ev, replyTo) => replyTo ! MetadataWriteSuccess(ev) }
        case Failure(regerts) =>
          log.error(regerts, "Failed to properly flush metadata to database")
          self ! DbWriteComplete
          putWithResponse foreach { case(ev, replyTo) => replyTo ! MetadataWriteFailure(regerts, ev) }
      }
      stay using NoData
    case Event(DbWriteComplete, curData) =>
      log.debug("Flush of metadata events complete")
      goto(WaitingToWrite) using curData
    // When receiving a put&respond message, add it to the current data so that when flushing metadata events, we have
    // enough information to be able to send an acknowledgement of success/failure of metadata event writes to the original requester.
    case Event(PutMetadataActionAndRespond(events, replyTo), curData) =>
      stay using curData.addData(PutMetadataActionAndRespond(events, replyTo))
  }

  onTransition {
    case WaitingToWrite -> WritingToDb => self ! FlushBatchToDb
  }
}

object WriteMetadataActor {
  def props(batchSize: Int, flushRate: FiniteDuration): Props = Props(new WriteMetadataActor(batchSize, flushRate)).withDispatcher(ServiceDispatcher)

  sealed trait WriteMetadataActorMessage
  case object CheckPendingWrites extends WriteMetadataActorMessage with MetadataServiceAction
  case object HasPendingWrites extends WriteMetadataActorMessage
  case object NoPendingWrites extends WriteMetadataActorMessage
}
