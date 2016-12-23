package cromwell.services.metadata.impl

import akka.actor.{ActorLogging, LoggingFSM, Props}
import cats.data.NonEmptyVector
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.SingletonServicesStore
import cromwell.services.metadata.MetadataEvent
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.impl.WriteMetadataActor.{WriteMetadataActorData, WriteMetadataActorState}

import scala.concurrent.duration._
import scala.util.{Failure, Success}

final case class WriteMetadataActor(batchRate: Int, flushRate: FiniteDuration)
  extends LoggingFSM[WriteMetadataActorState, WriteMetadataActorData] with ActorLogging with
  MetadataDatabaseAccess with SingletonServicesStore {
  import WriteMetadataActor._

  implicit val ec = context.dispatcher

  override def preStart(): Unit = {
    context.system.scheduler.schedule(0.seconds, flushRate, self, ScheduledFlushToDb)
    super.preStart()
  }

  startWith(WaitingToWrite, NoEvents)

  when(WaitingToWrite) {
    case Event(PutMetadataAction(events), curData) =>
      curData.addEvents(events) match {
        case data@HasEvents(e) if e.length > batchRate => goto(WritingToDb) using data
        case e => stay using e
      }
    case Event(ScheduledFlushToDb, curData) =>
      log.debug("Initiating periodic metadata flush to DB")
      goto(WritingToDb) using curData
  }

  when(WritingToDb) {
    case Event(ScheduledFlushToDb, curData) => stay using curData
    case Event(PutMetadataAction(events), curData) => stay using curData.addEvents(events)
    case Event(FlushBatchToDb, NoEvents) =>
      log.debug("Attempted metadata flush to DB but had nothing to write")
      goto(WaitingToWrite) using NoEvents
    case Event(FlushBatchToDb, HasEvents(e)) =>
      log.debug("Flushing {} metadata events to the DB", e.length)
      addMetadataEvents(e.toVector) onComplete {
        case Success(_) => self ! DbWriteComplete
        case Failure(regerts) =>
          log.error("Failed to properly flush metadata to database", regerts)
          self ! DbWriteComplete
      }

      stay using NoEvents
    case Event(DbWriteComplete, curData) =>
      log.debug("Flush of metadata events complete")
      goto(WaitingToWrite) using curData
  }

  onTransition {
    case WaitingToWrite -> WritingToDb => self ! FlushBatchToDb
  }
}

object WriteMetadataActor {
  def props(batchRate: Int, flushRate: FiniteDuration) = Props(new WriteMetadataActor(batchRate, flushRate)).withDispatcher(ServiceDispatcher)

  sealed trait WriteMetadataActorMessage
  case object DbWriteComplete extends WriteMetadataActorMessage
  case object FlushBatchToDb extends WriteMetadataActorMessage
  case object ScheduledFlushToDb extends WriteMetadataActorMessage

  sealed trait WriteMetadataActorState
  case object WaitingToWrite extends WriteMetadataActorState
  case object WritingToDb extends WriteMetadataActorState

  sealed trait WriteMetadataActorData {
    def addEvents(newEvents: Iterable[MetadataEvent]): WriteMetadataActorData = {
      NonEmptyVector.fromVector(newEvents.toVector) match {
        case Some(v) =>
          val newEvents = this match {
            case NoEvents => v
            case HasEvents(e) => e.concatNev(v)
          }
          HasEvents(newEvents)
        case None => this
      }
    }
  }

  case object NoEvents extends WriteMetadataActorData
  case class HasEvents(events: NonEmptyVector[MetadataEvent]) extends WriteMetadataActorData
}
