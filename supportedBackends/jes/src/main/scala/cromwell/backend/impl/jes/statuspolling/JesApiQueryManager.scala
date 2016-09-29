package cromwell.backend.impl.jes.statuspolling

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cromwell.backend.impl.jes.Run
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager._

/**
  * Currently, just holds a set of JES status poll requests until a PollingActor pulls the work.
  * TODO: Could eventually move all of the JES queries into a single work-pulling model.
  */
class JesApiQueryManager extends Actor with ActorLogging {

  private var workQueue: List[JesStatusPollQuery] = List.empty

  // TODO: Random magic number 5. Dynamically resize this?
  // TODO: Supervision strategy for pollers should be: stop that one, tell its requesters that it failed, and create a new one to replace it.
  private val pollers = (0 until 5) map { i => context.actorOf(JesPollingActor.props(self)) }
  pollers foreach { _ ! NoWorkToDo }

  override def receive = {
    case DoPoll(run) => workQueue :+= JesStatusPollQuery(sender, run)
    case RequestJesPollingWork =>
      log.debug(s"Request for JES Polling Work received (current queue size is ${workQueue.size})")
      handleJesPollingRequest(sender)
    case other => log.error(s"Unexpected message to JesPollingManager: $other")
  }

  /**
    * !! Function Impurity Warning: Modifies Actor Data !!
    */
  private def handleJesPollingRequest(workPullingJesPollingActor: ActorRef) = {
    val beheaded = beheadWorkQueue
    beheaded.workToDo match {
      case Some(work) =>
        log.debug(s"Sending work to JES poller.")
        workPullingJesPollingActor ! JesPollingWorkBatch(work)
      case None =>
        log.debug(s"No work for JES poller. Sad.")
        workPullingJesPollingActor ! NoWorkToDo
    }

    workQueue = beheaded.newWorkQueue
  }

  private final case class BeheadedWorkQueue(workToDo: Option[NonEmptyList[JesStatusPollQuery]], newWorkQueue: List[JesStatusPollQuery])
  private def beheadWorkQueue: BeheadedWorkQueue = {
    workQueue.grouped(100).toList match { // TODO: 100 is a random magic number: too big/too small/just right?
      case head :: tail => BeheadedWorkQueue(Option(NonEmptyList.fromListUnsafe(head)), tail.flatten) // NB: This IS safe: we know that head is non-empty due to the behavior of .grouped
      case Nil => BeheadedWorkQueue(None, List.empty)
    }
  }
}

object JesApiQueryManager {

  def props: Props = Props(new JesApiQueryManager)

  /**
    * Poll the job represented by the Run.
    */
  final case class DoPoll(run: Run)

  private[statuspolling] final case class JesStatusPollQuery(requester: ActorRef, run: Run)
  private[statuspolling] final case class JesPollingWorkBatch(workBatch: NonEmptyList[JesStatusPollQuery])
  private[statuspolling] case object NoWorkToDo

  private[statuspolling] case object RequestJesPollingWork
}
