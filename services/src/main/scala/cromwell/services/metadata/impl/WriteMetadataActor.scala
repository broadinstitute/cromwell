package cromwell.services.metadata.impl

import akka.actor.{ActorLogging, LoggingFSM, Props}
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.core.actor.BatchingDbWriter
import cromwell.core.actor.BatchingDbWriter._
import cromwell.services.SingletonServicesStore
import cromwell.services.metadata.MetadataEvent
import cromwell.services.metadata.MetadataService.{MetadataServiceAction, PutMetadataAction}

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
      val events = e.toVector.collect({ case e: MetadataEvent => e })
      addMetadataEvents(events) onComplete {
        case Success(_) => self ! DbWriteComplete
        case Failure(regerts) =>
          log.error(regerts, "Failed to properly flush metadata to database")
          self ! DbWriteComplete
      }

      stay using NoData
    case Event(DbWriteComplete, curData) =>
      log.debug("Flush of metadata events complete")
      goto(WaitingToWrite) using curData
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
