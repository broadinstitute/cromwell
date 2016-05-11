package cromwell.engine.workflow

import akka.actor.FSM.{CurrentState, SubscribeTransitionCallBack, Transition}
import akka.actor._
import akka.event.Logging
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.WorkflowId
import cromwell.engine._
import cromwell.engine.backend._
import cromwell.engine.db.DataAccess._
import cromwell.engine.workflow.WorkflowActor._
import cromwell.engine.workflow.WorkflowManagerActor.{AbortWorkflowCommand, _}
import cromwell.webservice.CromwellApiHandler._
import lenthall.config.ScalaConfig.EnhancedScalaConfig

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration._
import scala.concurrent.{Await, Future, Promise}
import scala.language.postfixOps

object WorkflowManagerActor {

  case class WorkflowIdToActorRef(workflowId: WorkflowId, workflowActor: ActorRef)

  sealed trait WorkflowManagerActorMessage
  /**
    * Commands
    */
  sealed trait WorkflowManagerActorCommand extends WorkflowManagerActorMessage

  case class SubmitWorkflowCommand(source: WorkflowSourceFiles) extends WorkflowManagerActorCommand
  case class AbortWorkflowCommand(id: WorkflowId) extends WorkflowManagerActorCommand
  case object AbortAllWorkflowsCommand extends WorkflowManagerActorCommand
  private[WorkflowManagerActor] final case class RestartWorkflowsCommand(workflows: Seq[OldStyleWorkflowDescriptor]) extends WorkflowManagerActorCommand

  /**
    * Responses
    */
  sealed trait WorkflowManagerActorResponse extends WorkflowManagerActorMessage

  def props(): Props = Props(new WorkflowManagerActor())

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
    def idFromActor(actor: ActorRef): Option[WorkflowId] =
      workflows.collectFirst { case (id, a) if a == actor => id }
    def withAddition(entry: WorkflowIdToActorRef): WorkflowManagerData =
      this.copy(workflows = workflows + (entry.workflowId -> entry.workflowActor))
    def without(id: WorkflowId): WorkflowManagerData =
      this.copy(workflows = workflows - id)
    def without(actor: ActorRef): WorkflowManagerData = {
      // If the ID was found in the lookup return a modified copy of the state data, otherwise just return
      // the same state data.
      idFromActor(actor) map without getOrElse this
    }
  }
}

class WorkflowManagerActor(config: Config)
  extends LoggingFSM[WorkflowManagerState, WorkflowManagerData] with CromwellActor {

  def this() = this(ConfigFactory.load)
  implicit val actorSystem = context.system

  private val RestartDelay: FiniteDuration = 200 milliseconds
  private val logger = Logging(context.system, this)
  private val tag = self.path.name

  private val donePromise = Promise[Unit]()

  private val conf = ConfigFactory.load

  override def preStart() {
    addShutdownHook()
    if (conf.getBoolean("system.workflow-restart")) { restartIncompleteWorkflows() }
  }

  private def addShutdownHook(): Unit = {
    // Only abort jobs on SIGINT if the config explicitly sets backend.abortJobsOnTerminate = true.
    val abortJobsOnTerminate =
      config.getConfig("backend").getBooleanOr("abortJobsOnTerminate", default = false)

    if (abortJobsOnTerminate) {
      sys.addShutdownHook {
        logger.info(s"$tag: Received shutdown signal. Aborting all running workflows...")
        self ! AbortAllWorkflowsCommand
        Await.ready(donePromise.future, Duration.Inf)
      }
    }
  }

  startWith(Running, WorkflowManagerData(workflows = Map.empty))

  when (Running) {
    /*
     Commands from clients
     */
    case Event(SubmitWorkflowCommand(source), stateData) =>
      val updatedEntry = submitWorkflow(source, replyTo = Option(sender), id = None)
      // Submit success brought in from the Cromwell API handler
      sender ! WorkflowManagerSubmitSuccess(updatedEntry.workflowId)
      stay() using stateData.withAddition(updatedEntry)
    case Event(AbortWorkflowCommand(id), stateData) =>
      val workflowActor = stateData.workflows.get(id)
      workflowActor match {
        case Some(actor) =>
          actor ! AbortWorkflowCommand
          stay()
        case None =>
          sender ! WorkflowManagerAbortFailure(id, new Exception(s"Couldn't abort $id because no workflow with that ID is in progress"))
          stay()
      }
    case Event(AbortAllWorkflowsCommand, data) if data.workflows.isEmpty =>
      goto(Done)
    case Event(AbortAllWorkflowsCommand, data) =>
      data.workflows.values.foreach { _ ! OldStyleWorkflowActor.AbortWorkflow }
      goto(Aborting)
    /*
     Internal commands
     */
    case Event(RestartWorkflowsCommand(w :: ws), stateData) =>
      val updatedEntry = restartWorkflow(w)
      context.system.scheduler.scheduleOnce(RestartDelay) { self ! RestartWorkflowsCommand(ws) }
      stay() using stateData.withAddition(updatedEntry)
    case Event(RestartWorkflowsCommand(Nil), _) =>
      // No more workflows need restarting.
      stay()
    /*
     Responses from services
     */
    case Event(WorkflowSucceededResponse(workflowId), data) =>
      log.info(s"Workflow $workflowId succeeded!")
      stay()
    case Event(WorkflowFailedResponse(workflowId, inState, reasons), data) =>
      log.error(s"Workflow $workflowId failed (during $inState): ${reasons.mkString("\n")}")
      stay()
    /*
     Watched transitions
     */
    case Event(Transition(workflowActor, _, toState: WorkflowActorTerminalState), data) =>
      log.info(workflowActor.path.name + "has gone terminal")
      stay using data.without(workflowActor)
  }

  when (Aborting) {
    case Event(Transition(workflowActor, _, toState: WorkflowState), data) if toState.isTerminal =>
      // Remove this terminal actor from the workflowStore and log a progress message.
      val updatedData = data.without(workflowActor)
      logger.info(s"$tag: Waiting for all workflows to abort (${updatedData.workflows.size} remaining).")
      // If there are no more workflows to abort we're done, otherwise just stay in the current state.
      val resultAction = if (updatedData.workflows.isEmpty) goto(Done) else stay()
      // Whatever the result action, use the updated data:
      resultAction using updatedData
  }

  when (Done) { FSM.NullFunction }

  whenUnhandled {
    // Uninteresting transition and current state notifications.
    case Event((Transition(_, _, _) | CurrentState(_, _)), _) => stay()
    // Anything else certainly IS interesting:
    case Event(unhandled, data) =>
      log.warning(s"$tag: Unhandled message: $unhandled")
      stay()
  }

  onTransition {
    case _ -> Done =>
      logger.info(s"$tag: All workflows finished. Shutting down.")
      donePromise.trySuccess(())
    case fromState -> toState =>
      logger.info(s"$tag transitioning from $fromState to $toState")
  }


  /** Submit the workflow and return an updated copy of the state data reflecting the addition of a
    * Workflow ID -> WorkflowActorRef entry.
    */
  private def submitWorkflow(source: WorkflowSourceFiles, replyTo: Option[ActorRef], id: Option[WorkflowId]): WorkflowIdToActorRef = {
    val workflowId: WorkflowId = id.getOrElse(WorkflowId.randomId())
    logger.info(s"$tag submitWorkflow input id = $id, effective id = $workflowId")
    val isRestart = id.isDefined

    val startMode = if (isRestart) RestartExistingWorkflow else StartNewWorkflow
    val wfActor = context.actorOf(WorkflowActor.props(workflowId, startMode, source, config), name = s"WorkflowActor-$workflowId")
    replyTo.foreach { _ ! WorkflowManagerSubmitSuccess(id = workflowId) }
    wfActor ! SubscribeTransitionCallBack(self)
    wfActor ! StartWorkflowCommand
    logger.info(s"Successfuly started ${wfActor.path} for Workflow $workflowId")
    WorkflowIdToActorRef(workflowId, wfActor)
  }

  private def restartWorkflow(restartableWorkflow: OldStyleWorkflowDescriptor): WorkflowIdToActorRef = {
    logger.info("Invoking restartableWorkflow on " + restartableWorkflow.id.shortString)
    submitWorkflow(restartableWorkflow.sourceFiles, replyTo = None, Option(restartableWorkflow.id))
  }

  private def restartIncompleteWorkflows(): Unit = {
    def logRestarts(restartableWorkflows: Traversable[OldStyleWorkflowDescriptor]): Unit = {
      val num = restartableWorkflows.size
      val displayNum = if (num == 0) "no" else num.toString
      val plural = if (num == 1) "" else "s"

      logger.info(s"$tag Found $displayNum workflow$plural to restart.")

      if (num > 0) {
        val ids = restartableWorkflows.map { _.id.toString }.toSeq.sorted
        logger.info(s"$tag Restarting workflow ID$plural: " + ids.mkString(", "))
      }
    }

    // TODO: This floating Future is awkward in an actor method. I would like to move all of this out into a WorkflowRestarterActor
    val result = for {
      restartableWorkflowExecutionAndAuxes <- globalDataAccess.getWorkflowExecutionAndAuxByState(Seq(WorkflowSubmitted, WorkflowRunning))
      restartableWorkflows <- Future.sequence(restartableWorkflowExecutionAndAuxes map workflowDescriptorFromExecutionAndAux)
      _ = logRestarts(restartableWorkflows)
      _ = self ! RestartWorkflowsCommand(restartableWorkflows.toSeq)
    } yield ()

    result recover {
      case e: Throwable => logger.error(e, e.getMessage)
    }
  }
}
