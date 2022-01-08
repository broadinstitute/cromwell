package cromwell.core.actor

import java.time.OffsetDateTime

import akka.actor.{ActorRef, FSM, Timers}
import akka.dispatch.ControlMessage
import cats.data.NonEmptyVector
import common.collections.WeightedQueue
import cromwell.core.actor.BatchActor._
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.slf4j.LoggerFactory

import scala.concurrent.Future
import scala.concurrent.duration.{Duration, FiniteDuration}
import scala.util.{Failure, Success}


/** A collection of state, data, and message types to support BatchActor. */
object BatchActor {
  type BatchData[C] = WeightedQueue[C, Int]

  val logger = LoggerFactory.getLogger("BatchActor")

  /** The states for BatchActor. */
  sealed trait BatchActorState
  case object WaitingToProcess extends BatchActorState
  case object Processing extends BatchActorState

  /**
    * The control messages for BatchActor. Because they extend ControlMessage, in combination with a
    * control aware mailbox those messages will be processed in priority.
    */
  sealed trait BatchActorControlMessage extends ControlMessage
  case object ProcessingComplete extends BatchActorControlMessage
  case object ScheduledFlushKey
  case object ScheduledProcessAction extends BatchActorControlMessage

  case class CommandAndReplyTo[C](command: C, replyTo: ActorRef)

  case class QueueWeight(weight: Int)
}

/**
  * Abstract FSM that provides queueing and batching behavior based on a batch size and a flush rate.
  * It is backed by a WeightedQueue which makes it possible to decouple the number of messages received from
  * the effective "weight" of the queue.
  */
abstract class BatchActor[C](val flushRate: FiniteDuration,
                             val batchSize: Int) extends FSM[BatchActorState, BatchData[C]] with Timers {
  private var shuttingDown: Boolean = false

  implicit val ec = context.dispatcher
  private val name = self.path.name

  /**
    * Override to false to prevent this actor from logging its configuration on startup.
    */
  protected def logOnStartUp: Boolean = true

  /**
    * Override to true if this actor is going to be placed behind a router.
    */
  protected def routed: Boolean = false

  override def preStart(): Unit = {
    if (logOnStartUp) log.info("{} configured to flush with batch size {} and process rate {}.", name, batchSize, flushRate)
    if (flushRate != Duration.Zero) {
      timers.startPeriodicTimer(ScheduledFlushKey, ScheduledProcessAction, flushRate)
    }
    super.preStart()
  }

  startWith(WaitingToProcess, WeightedQueue.empty[C, Int](weightFunction))

  def commandToData(snd: ActorRef): PartialFunction[Any, C]

  // If work has arrived 'recently' - ie within this duration - ignore any periodic 'ScheduledProcessAction' messages
  val recentArrivalThreshold: Option[FiniteDuration] = None
  var mostRecentArrival: Option[OffsetDateTime] = None
  def suitableIntervalSinceLastArrival(): Boolean = (for {
    threshold <- recentArrivalThreshold
    mostRecent <- mostRecentArrival
    now = OffsetDateTime.now
  } yield mostRecent.plusNanos(threshold.toNanos).isBefore(now)).getOrElse(true)

  when(WaitingToProcess) {
    // On a regular event, only process if the batch size has been reached.
    case Event(command, data) if commandToData(sender).isDefinedAt(command) =>
      recentArrivalThreshold foreach { _ => mostRecentArrival = Option(OffsetDateTime.now()) }
      processIfBatchSizeReached(data.enqueue(commandToData(sender)(command)))
    // On a scheduled process, always process
    case Event(ScheduledProcessAction, data) =>
      if (suitableIntervalSinceLastArrival()) {
        gossip(QueueWeight(data.weight))
        processHead(data)
      } else {
        stay()
      }
    case Event(ShutdownCommand, data) =>
      logger.info(s"{} Shutting down: ${data.weight} queued messages to process", self.path.name)
      shuttingDown = true
      processHead(data)
  }

  when(Processing) {
    // Already processing, enqueue the command
    case Event(command, data) if commandToData(sender).isDefinedAt(command) =>
      stay() using data.enqueue(commandToData(sender)(command))
    // Already processing, can only do one at a time
    case Event(ScheduledProcessAction, data) => 
      gossip(QueueWeight(data.weight))
      stay()
    // Process is complete and we're shutting down so process even if we're under the batch size.
    case Event(ProcessingComplete, data) if shuttingDown =>
      logger.info(s"{} Shutting down: processing ${data.weight} queued messages", self.path.name)
      processHead(data)
    // Processing is complete, re-process only if needed    
    case Event(ProcessingComplete, data) if !shuttingDown =>
      processIfBatchSizeReached(data)
    case Event(ShutdownCommand, _) =>
      shuttingDown = true
      stay()
  }

  /**
    * Given a command, return its weight
    */
  protected def weightFunction(command: C): Int

  /**
    * Process the data asynchronously
    * @return the number of elements processed
    */
  protected def process(data: NonEmptyVector[C]): Future[Int]

  private def processIfBatchSizeReached(data: BatchData[C]) = {
    if (data.weight >= batchSize) processHead(data)
    else goto(WaitingToProcess) using data
  }

  private def processHead(data: BatchData[C]) = if (data.innerQueue.nonEmpty) {
    val (head, newQueue) = data.behead(batchSize)

    def processNonEmptyHead(nev: NonEmptyVector[C]) = process(nev) onComplete {
      case Success(_) =>
        self ! ProcessingComplete
      case Failure(regerts) =>
        log.error(regerts, "{} Failed to properly process data", name)
        self ! ProcessingComplete
    }

    head.headOption match {
      case Some(headOfHead) => 
        processNonEmptyHead(NonEmptyVector(headOfHead, head.tail))
        goto(Processing) using newQueue
      case None =>
        goto(WaitingToProcess) using data
    }

    goto(Processing) using newQueue
  } else
  // This goto is important, even if we're already in WaitingToProcess we want to trigger the onTransition below
  // to check if it's time to shutdown
    goto(WaitingToProcess) using data

  onTransition {
    case _ -> WaitingToProcess if shuttingDown && nextStateData.innerQueue.isEmpty =>
      timers.cancelAll()
      context stop self
  }
}
