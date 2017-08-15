package cromwell.core.actor

import akka.actor.{ActorRef, Cancellable, FSM}
import cats.data.NonEmptyVector
import cromwell.core.actor.BatchingDbWriter._
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.slf4j.LoggerFactory

import scala.util.{Failure, Success, Try}


/** A collection of state, data, and message types to support batched database writes. */
object BatchingDbWriter {
  val logger = LoggerFactory.getLogger("BatchingDbWriteActor")

  /** Data for batched database writes. */
  sealed trait BatchingDbWriterData {
    def addData[D](datum: D): BatchingDbWriterData = addData(Vector(datum))
    def addData[D](data: Iterable[D]): BatchingDbWriterData = {
      Try(NonEmptyVector.fromVector(data.toVector)) match {
        case Success(Some(v)) =>
          val newEvents = this match {
            case NoData => v
            case HasData(e) => e.concatNev(v)
          }
          HasData(newEvents)
        case Success(None) => this
        case Failure(f) =>
          val dataSample = data.take(3).mkString(", ") + (if (data.size > 3) ", ..." else "")
          logger.error(s"Failed processing batched data: $dataSample. Data will be dropped and not be sent to the database.", f)
          this
      }
    }

    def length: Int = this match {
      case NoData => 0
      case HasData(e) => e.length
    }
  }

  case object NoData extends BatchingDbWriterData
  case class HasData[E](events: NonEmptyVector[E]) extends BatchingDbWriterData

  /** The states for batched database writes. */
  sealed trait BatchingDbWriterState
  case object WaitingToWrite extends BatchingDbWriterState
  case object WritingToDb extends BatchingDbWriterState

  /** The message types for batched database writes. */
  sealed trait BatchingDbWriterMessage
  case object DbWriteComplete extends BatchingDbWriterMessage
  case object FlushBatchToDb extends BatchingDbWriterMessage
  case object ScheduledFlushToDb extends BatchingDbWriterMessage

  case class CommandAndReplyTo[C](command: C, replyTo: ActorRef)
}

/**
  * Trait that contains some common batch-related and graceful shutdown logic.
  * Be careful NOT to add a custom whenUnhandled state function when mixing in this trait as it will override the 
  * graceful shutdown handling logic.
  * 
  * Note that there is more common logic that could be abstracted here.
  */
trait BatchingDbWriterActor { this: FSM[BatchingDbWriterState, BatchingDbWriterData] =>
  import scala.concurrent.duration._
  
  private var shuttingDown: Boolean = false
  
  def isShuttingDown: Boolean = shuttingDown
  def dbFlushRate: FiniteDuration
  def dbBatchSize: Int

  var periodicFlush: Option[Cancellable] = None

  override def preStart(): Unit = {
    log.info("{} configured to write to the database with batch size {} and flush rate {}.", self.path.name, dbBatchSize, dbFlushRate)
    periodicFlush = Option(context.system.scheduler.schedule(0.seconds, dbFlushRate, self, ScheduledFlushToDb)(context.dispatcher))
  }

  /**
    * WhenUnhandled state function that handles reception of ShutdownCommand and acts appropriately
    */
  private val whenUnhandledFunction: StateFunction = {
    case Event(ShutdownCommand, NoData) if stateName == WaitingToWrite =>
      periodicFlush foreach { _.cancel() }
      context stop self
      stay()
    case Event(ShutdownCommand, _) if stateName == WaitingToWrite =>
      logger.info("{} flushing database writes...", self.path.name)
      shuttingDown = true
      // transitioning to WritingToDb triggers a FlushBatchToDb to be sent to self
      goto(WritingToDb)
    case Event(ShutdownCommand, _) if stateName == WritingToDb =>
      logger.info("{} waiting for database writes to be flushed...", self.path.name)
      shuttingDown = true
      stay()
  }

  whenUnhandled(whenUnhandledFunction)

  onTransition {
    case WaitingToWrite -> WritingToDb => self ! FlushBatchToDb
    // When transitioning back to WaitingToWrite, if there's no data left to process, and we're trying to shutdown, then stop
    case _ -> WaitingToWrite if shuttingDown && nextStateData == NoData =>
      periodicFlush foreach { _.cancel() }
      context stop self
  }
}
