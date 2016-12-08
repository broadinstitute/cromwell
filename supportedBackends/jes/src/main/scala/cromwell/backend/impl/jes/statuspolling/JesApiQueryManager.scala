package cromwell.backend.impl.jes.statuspolling

import akka.actor.{Actor, ActorLogging, ActorRef, Props, SupervisorStrategy, Terminated}
import cats.data.NonEmptyList
import cromwell.backend.impl.jes.Run
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager._
import cromwell.core.Dispatcher.BackendDispatcher

import scala.collection.immutable.Queue

/**
  * Currently, just holds a set of JES status poll requests until a PollingActor pulls the work.
  * TODO: Could eventually move all of the JES queries into a single work-pulling model.
  */
class JesApiQueryManager(val qps: Int) extends Actor with ActorLogging {

  // workQueue is protected for the unit tests, not intended to be generally overridden
  protected[statuspolling] var workQueue: Queue[JesStatusPollQuery] = Queue.empty
  private var workInProgress: Map[ActorRef, JesPollingWorkBatch] = Map.empty

  // If the statusPoller dies, we want to stop it and handle the termination ourselves.
  override val supervisorStrategy = SupervisorStrategy.stoppingStrategy
  private def statusPollerProps = JesPollingActor.props(self, qps)

  // statusPoller is protected for the unit tests, not intended to be generally overridden
  protected[statuspolling] var statusPoller: ActorRef = _

  resetStatusPoller()

  override def receive = {
    case DoPoll(run) => workQueue :+= JesStatusPollQuery(sender, run)
    case RequestJesPollingWork(maxBatchSize) =>
      log.debug(s"Request for JES Polling Work received (max batch: $maxBatchSize, current queue size is ${workQueue.size})")
      handleJesPollingRequest(sender, maxBatchSize)
    case Terminated(actorRef) => handleTerminated(actorRef)
    case other => log.error(s"Unexpected message to JesPollingManager: $other")
  }

  /**
    * !! Function Impurity Warning: Modifies Actor Data !!
    */
  private def handleJesPollingRequest(workPullingJesPollingActor: ActorRef, maxBatchSize: Int) = {
    workInProgress -= workPullingJesPollingActor
    val beheaded = beheadWorkQueue(maxBatchSize)
    beheaded.workToDo match {
      case Some(work) =>
        sendWork(workPullingJesPollingActor, JesPollingWorkBatch(work))
      case None =>
        log.debug(s"No work for JES poller. Sad.")
        workPullingJesPollingActor ! NoWorkToDo
    }

    workQueue = beheaded.newWorkQueue
  }

  private def sendWork(workPullingJesPollingActor: ActorRef, work: JesPollingWorkBatch) = {
    log.debug(s"Sending work to JES poller.")
    workPullingJesPollingActor ! work
    workInProgress += (workPullingJesPollingActor -> work)
  }

  private final case class BeheadedWorkQueue(workToDo: Option[NonEmptyList[JesStatusPollQuery]], newWorkQueue: Queue[JesStatusPollQuery])
  private def beheadWorkQueue(maxBatchSize: Int): BeheadedWorkQueue = {

    val head = workQueue.take(maxBatchSize).toList
    val tail = workQueue.drop(maxBatchSize)

    head match {
      case h :: t => BeheadedWorkQueue(Option(NonEmptyList(h, t)), tail)
      case Nil => BeheadedWorkQueue(None, Queue.empty)
    }
  }

  private def handleTerminated(terminee: ActorRef) = {
    // Currently we can assume this is a polling actor. Might change in a future update:
    workInProgress.get(terminee) match {
      case Some(work) =>
        // Most likely due to an unexpected HTTP error, push the work back on the queue and keep going
        log.info(s"The JES polling actor $terminee unexpectedly terminated while conducting ${work.workBatch.tail.size + 1} polls. Making a new one...")
        workInProgress -= terminee
        workQueue = workQueue ++ work.workBatch.toList
      case None =>
        // It managed to die while doing absolutely nothing...!?
        // Maybe it deserves an entry in https://en.wikipedia.org/wiki/List_of_unusual_deaths
        // Oh well, in the mean time don't do anything, just start a new one
        log.error(s"The JES polling actor $terminee managed to unexpectedly terminate whilst doing absolutely nothing. This is probably a programming error. Making a new one...")
    }

    resetStatusPoller()
  }

  private def resetStatusPoller() = {
    statusPoller = makeStatusPoller()
    context.watch(statusPoller)
    log.info(s"watching $statusPoller")
  }

  private[statuspolling] def makeStatusPoller(): ActorRef = context.actorOf(statusPollerProps)
}

object JesApiQueryManager {

  def props(qps: Int): Props = Props(new JesApiQueryManager(qps)).withDispatcher(BackendDispatcher)

  /**
    * Poll the job represented by the Run.
    */
  final case class DoPoll(run: Run)

  private[statuspolling] final case class JesStatusPollQuery(requester: ActorRef, run: Run)
  private[statuspolling] final case class JesPollingWorkBatch(workBatch: NonEmptyList[JesStatusPollQuery])
  private[statuspolling] case object NoWorkToDo

  private[statuspolling] final case class RequestJesPollingWork(maxBatchSize: Int)
}
