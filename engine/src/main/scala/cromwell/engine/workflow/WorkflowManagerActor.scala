package cromwell.engine.workflow

import akka.actor.FSM.{CurrentState, SubscribeTransitionCallBack, Transition}
import akka.actor._
import akka.event.Logging
import cats.data.NonEmptyList
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.{WorkflowAborted, WorkflowId}
import cromwell.engine.backend.BackendSingletonCollection
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.WorkflowManagerActor._
import cromwell.engine.workflow.workflowstore.{WorkflowStoreActor, WorkflowStoreState}
import cromwell.jobstore.JobStoreActor.{JobStoreWriteFailure, JobStoreWriteSuccess, RegisterWorkflowCompleted}
import cromwell.services.metadata.MetadataService._
import cromwell.webservice.EngineStatsActor
import net.ceedubs.ficus.Ficus._
import org.apache.commons.lang3.exception.ExceptionUtils

import scala.concurrent.duration._

object WorkflowManagerActor {
  val DefaultMaxWorkflowsToRun = 5000
  val DefaultMaxWorkflowsToLaunch = 50
  val DefaultNewWorkflowPollRate = 20

  case class WorkflowIdToActorRef(workflowId: WorkflowId, workflowActor: ActorRef)
  class WorkflowNotFoundException(s: String) extends Exception(s)

  sealed trait WorkflowManagerActorMessage
  /**
    * Commands
    */
  sealed trait WorkflowManagerActorCommand extends WorkflowManagerActorMessage
  case object RetrieveNewWorkflows extends WorkflowManagerActorCommand
  final case class AbortWorkflowCommand(id: WorkflowId, replyTo: ActorRef) extends WorkflowManagerActorCommand
  case object AbortAllWorkflowsCommand extends WorkflowManagerActorCommand
  final case class SubscribeToWorkflowCommand(id: WorkflowId) extends WorkflowManagerActorCommand
  case object EngineStatsCommand extends WorkflowManagerActorCommand

  def props(workflowStore: ActorRef,
            serviceRegistryActor: ActorRef,
            workflowLogCopyRouter: ActorRef,
            jobStoreActor: ActorRef,
            callCacheReadActor: ActorRef,
            jobTokenDispenserActor: ActorRef,
            backendSingletonCollection: BackendSingletonCollection): Props = {
    Props(new WorkflowManagerActor(
      workflowStore, serviceRegistryActor, workflowLogCopyRouter, jobStoreActor, callCacheReadActor, jobTokenDispenserActor, backendSingletonCollection)
    ).withDispatcher(EngineDispatcher)
  }

  /**
    * States
    */
  sealed trait WorkflowManagerState
  case object Running extends WorkflowManagerState
  case object Aborting extends WorkflowManagerState
  case object Done extends WorkflowManagerState

  /**
    * Data
    */
  case class WorkflowManagerData(workflows: Map[WorkflowId, ActorRef]) {
    def idFromActor(actor: ActorRef): Option[WorkflowId] = workflows.collectFirst { case (id, a) if a == actor => id }

    def withAddition(entries: NonEmptyList[WorkflowIdToActorRef]): WorkflowManagerData = {
      val entryTuples = entries map { e => e.workflowId -> e.workflowActor }
      this.copy(workflows = workflows ++ entryTuples.toList)
    }

    def without(id: WorkflowId): WorkflowManagerData = this.copy(workflows = workflows - id)
    def without(actor: ActorRef): WorkflowManagerData = {
      // If the ID was found in the lookup return a modified copy of the state data, otherwise just return
      // the same state data.
      idFromActor(actor) map without getOrElse this
    }
  }
}

class WorkflowManagerActor(config: Config,
                           val workflowStore: ActorRef,
                           val serviceRegistryActor: ActorRef,
                           val workflowLogCopyRouter: ActorRef,
                           val jobStoreActor: ActorRef,
                           val callCacheReadActor: ActorRef,
                           val jobTokenDispenserActor: ActorRef,
                           val backendSingletonCollection: BackendSingletonCollection)
  extends LoggingFSM[WorkflowManagerState, WorkflowManagerData] {

  def this(workflowStore: ActorRef,
           serviceRegistryActor: ActorRef,
           workflowLogCopyRouter: ActorRef,
           jobStoreActor: ActorRef,
           callCacheReadActor: ActorRef,
           jobTokenDispenserActor: ActorRef,
           backendSingletonCollection: BackendSingletonCollection) = this(
    ConfigFactory.load, workflowStore, serviceRegistryActor, workflowLogCopyRouter, jobStoreActor, callCacheReadActor, jobTokenDispenserActor, backendSingletonCollection)

  private val maxWorkflowsRunning = config.getConfig("system").as[Option[Int]]("max-concurrent-workflows").getOrElse(DefaultMaxWorkflowsToRun)
  private val maxWorkflowsToLaunch = config.getConfig("system").as[Option[Int]]("max-workflow-launch-count").getOrElse(DefaultMaxWorkflowsToLaunch)
  private val newWorkflowPollRate = config.getConfig("system").as[Option[Int]]("new-workflow-poll-rate").getOrElse(DefaultNewWorkflowPollRate).seconds

  private val logger = Logging(context.system, this)
  private val tag = self.path.name

  private var abortingWorkflowToReplyTo = Map.empty[WorkflowId, ActorRef]

  override def preStart(): Unit = {
    addShutdownHook()
    // Starts the workflow polling cycle
    self ! RetrieveNewWorkflows
  }

  private def addShutdownHook() = {
    // Only abort jobs on SIGINT if the config explicitly sets system.abort-jobs-on-terminate = true.
    val abortJobsOnTerminate =
      config.getConfig("system").as[Option[Boolean]]("abort-jobs-on-terminate").getOrElse(false)

    if (abortJobsOnTerminate) {
      sys.addShutdownHook {
        logger.info(s"$tag: Received shutdown signal.")
        self ! AbortAllWorkflowsCommand
        while (stateData != null && stateData.workflows.nonEmpty) {
          log.info(s"Waiting for ${stateData.workflows.size} workflows to abort...")
          Thread.sleep(1000)
        }
      }
    }
  }

  startWith(Running, WorkflowManagerData(workflows = Map.empty))

  when (Running) {
    /*
     Commands from clients
     */
    case Event(RetrieveNewWorkflows, stateData) =>
      /*
        Cap the total number of workflows in flight, but also make sure we don't pull too many in at once.
        Determine the number of available workflow slots and request the smaller of that number of maxWorkflowsToLaunch.
       */
      val maxNewWorkflows = maxWorkflowsToLaunch min (maxWorkflowsRunning - stateData.workflows.size)
      workflowStore ! WorkflowStoreActor.FetchRunnableWorkflows(maxNewWorkflows)
      stay()
    case Event(WorkflowStoreActor.NoNewWorkflowsToStart, stateData) =>
      log.debug("WorkflowStore provided no new workflows to start")
      scheduleNextNewWorkflowPoll()
      stay()
    case Event(WorkflowStoreActor.NewWorkflowsToStart(newWorkflows), stateData) =>
      val newSubmissions = newWorkflows map submitWorkflow
      log.info("Retrieved {} workflows from the WorkflowStoreActor", newSubmissions.toList.size)
      scheduleNextNewWorkflowPoll()
      stay() using stateData.withAddition(newSubmissions)
    case Event(SubscribeToWorkflowCommand(id), data) =>
      data.workflows.get(id) foreach {_ ! SubscribeTransitionCallBack(sender())}
      stay()
    case Event(WorkflowManagerActor.AbortWorkflowCommand(id, replyTo), stateData) =>
      val workflowActor = stateData.workflows.get(id)
      workflowActor match {
        case Some(actor) =>
          actor ! WorkflowActor.AbortWorkflowCommand
          abortingWorkflowToReplyTo += id -> replyTo
          // Wait until the workflow transitions to the aborted state to respond with `WorkflowAborted`.
          stay()
        case None =>
          // All cool, if we got this far the workflow ID was found in the workflow store so this workflow must have never
          // made it to the workflow manager.
          replyTo ! WorkflowStoreActor.WorkflowAborted(id)
          stay()
      }
    case Event(AbortAllWorkflowsCommand, data) if data.workflows.isEmpty =>
      goto(Done)
    case Event(AbortAllWorkflowsCommand, data) =>
      log.info(s"$tag Aborting all workflows")
      data.workflows.values.foreach { _ ! WorkflowActor.AbortWorkflowCommand }
      goto(Aborting)
    /*
     Responses from services
     */
    case Event(WorkflowFailedResponse(workflowId, inState, reasons), data) =>
      log.error(s"$tag Workflow $workflowId failed (during $inState): ${expandFailureReasons(reasons)}")
      stay()
    /*
     Watched transitions
     */
    case Event(Transition(workflowActor, _, toState: WorkflowActorTerminalState), data) =>
      log.info(s"$tag ${workflowActor.path.name} is in a terminal state: $toState")
      // This silently fails if idFromActor is None, but data.without call right below will as well
      data.idFromActor(workflowActor) foreach { workflowId =>
        jobStoreActor ! RegisterWorkflowCompleted(workflowId)
        if (toState.workflowState == WorkflowAborted) {
          val replyTo = abortingWorkflowToReplyTo(workflowId)
          replyTo ! WorkflowStoreActor.WorkflowAborted(workflowId)
          abortingWorkflowToReplyTo -= workflowId
        } else {
          workflowStore ! WorkflowStoreActor.RemoveWorkflow(workflowId)
        }
      }
      stay using data.without(workflowActor)
  }

  when (Aborting) {
    case Event(Transition(workflowActor, _, toState: WorkflowActorState), data) if toState.terminal =>
      // Remove this terminal actor from the workflowStore and log a progress message.
      val updatedData = data.without(workflowActor)
      // If there are no more workflows to abort we're done, otherwise just stay in the current state.
      val resultAction = if (updatedData.workflows.isEmpty) {
        logger.info(s"$tag All workflows are aborted")
        goto(Done)
      } else {
        logger.info(s"$tag Waiting for all workflows to abort (${updatedData.workflows.size} remaining).")
        stay()
      }
      // Whatever the result action, use the updated data:
      resultAction using updatedData
    case Event(_, _) => stay()
  }

  when (Done) { FSM.NullFunction }

  whenUnhandled {
    case Event(MetadataPutFailed(action, error), _) =>
      log.warning(s"$tag Put failed for Metadata action $action : ${error.getMessage}")
      stay()
    case Event(MetadataPutAcknowledgement(_), _) => stay()
    // Uninteresting transition and current state notifications.
    case Event((Transition(_, _, _) | CurrentState(_, _)), _) => stay()
    case Event(JobStoreWriteSuccess(_), _) => stay() // Snoozefest
    case Event(JobStoreWriteFailure(t), _) =>
      log.error("Error writing to JobStore from WorkflowManagerActor: {}", t)
      // This is minorly bad. The JobStore would never rid itself of this workflow's entries. Unlikely to be a big deal
      stay()
    case Event(EngineStatsCommand, data) =>
      val sndr = sender()
      context.actorOf(EngineStatsActor.props(data.workflows.values.toList, sndr), s"EngineStatsActor-${sndr.hashCode()}")
      stay()
    // Anything else certainly IS interesting:
    case Event(unhandled, data) =>
      log.warning(s"$tag Unhandled message: $unhandled")
      stay()
  }

  onTransition {
    case _ -> Done =>
      logger.info(s"$tag All workflows finished")
      ()
    case fromState -> toState =>
      logger.debug(s"$tag transitioning from $fromState to $toState")
  }

  /**
    * Submit the workflow and return an updated copy of the state data reflecting the addition of a
    * Workflow ID -> WorkflowActorRef entry.
    */
  private def submitWorkflow(workflow: workflowstore.WorkflowToStart): WorkflowIdToActorRef = {
    val workflowId = workflow.id

    val startMode = if (workflow.state == WorkflowStoreState.Restartable) {
      logger.info(s"$tag Restarting workflow UUID($workflowId)")
      RestartExistingWorkflow
    } else {
      logger.info(s"$tag Starting workflow UUID($workflowId)")
      StartNewWorkflow
    }

    val wfProps = WorkflowActor.props(workflowId, startMode, workflow.sources, config, serviceRegistryActor,
      workflowLogCopyRouter, jobStoreActor, callCacheReadActor, jobTokenDispenserActor, backendSingletonCollection)
    val wfActor = context.actorOf(wfProps, name = s"WorkflowActor-$workflowId")

    wfActor ! SubscribeTransitionCallBack(self)
    wfActor ! StartWorkflowCommand
    logger.info(s"$tag Successfully started ${wfActor.path.name}")
    WorkflowIdToActorRef(workflowId, wfActor)
  }

  private def scheduleNextNewWorkflowPoll() = {
    context.system.scheduler.scheduleOnce(newWorkflowPollRate, self, RetrieveNewWorkflows)(context.dispatcher)
  }

  private def expandFailureReasons(reasons: Seq[Throwable]) = {
    reasons map { reason =>
      reason.getMessage + "\n" + ExceptionUtils.getStackTrace(reason)
    } mkString "\n"
  }
}
