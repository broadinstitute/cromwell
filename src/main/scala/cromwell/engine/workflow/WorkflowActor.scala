package cromwell.engine.workflow

import akka.actor.{FSM, LoggingFSM, Props}
import akka.event.Logging
import akka.pattern.pipe
import cromwell.binding._
import cromwell.binding.values.{WdlObject, WdlValue}
import cromwell.engine._
import cromwell.engine.ExecutionStatus._
import cromwell.engine.backend.Backend
import cromwell.engine.db.DataAccess.WorkflowInfo
import cromwell.engine.db.{CallStatus, DataAccess}
import cromwell.engine.workflow.WorkflowActor._

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration._
import scala.concurrent.{Await, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object WorkflowActor {
  sealed trait WorkflowActorMessage
  case object Start extends WorkflowActorMessage
  case object Restart extends WorkflowActorMessage
  case object Complete extends WorkflowActorMessage
  case object GetFailureMessage extends WorkflowActorMessage
  case object GetOutputs extends WorkflowActorMessage
  case class CallStarted(call: Call) extends WorkflowActorMessage
  case class CallCompleted(call: Call, callOutputs: CallOutputs) extends WorkflowActorMessage
  case class CallFailed(call: Call, failure: String) extends WorkflowActorMessage
  case class RunnableCalls(calls: Iterable[Call]) extends WorkflowActorMessage

  def props(descriptor: WorkflowDescriptor, backend: Backend, dataAccess: DataAccess): Props = {
    Props(WorkflowActor(descriptor, backend, dataAccess))
  }

  sealed trait WorkflowFailure
  case object NoFailureMessage extends WorkflowFailure
  case class FailureMessage(msg: String) extends WorkflowFailure with WorkflowActorMessage

  val DatabaseTimeout = 5 seconds

  type ExecutionStore = Map[Call, ExecutionStatus.Value]
}

case class WorkflowActor(workflow: WorkflowDescriptor,
                         backend: Backend,
                         dataAccess: DataAccess)
  extends LoggingFSM[WorkflowState, WorkflowFailure] with CromwellActor {

  private var executionStore: Map[Call, ExecutionStatus.Value] = _

  val tag: String = s"WorkflowActor [UUID(${workflow.shortId})]"
  override val log = Logging(context.system, classOf[WorkflowActor])

  startWith(WorkflowSubmitted, NoFailureMessage)

  def initWorkflow(initialization: Future[Unit] = Future.successful(())): ExecutionStore = {
    val futureStore = for {
      _ <- initialization
      store <- createStore
    } yield store
    Await.result(futureStore, DatabaseTimeout)
  }

  /**
   * Attempt to start all runnable calls.  If successful, go to the `WorkflowRunning` state if not already in that
   * state, otherwise remain in `WorkflowRunning`.  If not successful, go to `WorkflowFailed`.
   */
  private def startRunnableCalls() = {
      tryStartingRunnableCalls() match {
        case Success(_) =>
          goto(WorkflowRunning)
        case Failure(e) =>
          log.error(e, e.getMessage)
          goto(WorkflowFailed)
      }
  }

  when(WorkflowSubmitted) {
    case Event(Restart, NoFailureMessage) =>
      executionStore = initWorkflow()
      startRunnableCalls()
    case Event(Start, NoFailureMessage) =>
      log.info(s"$tag Start message received")
      executionStore = initWorkflow(createWorkflow())
      startRunnableCalls()
  }

  when(WorkflowRunning) {
    case Event(CallStarted(call), NoFailureMessage) =>
      persistStatus(call, ExecutionStatus.Running)
      stay()
    case Event(CallCompleted(call, outputs), NoFailureMessage) =>
      awaitCallComplete(call, outputs) match {
        case Success(_) =>
          if (isWorkflowDone) goto(WorkflowSucceeded) else startRunnableCalls()
        case Failure(e) =>
          log.error(e, e.getMessage)
          goto(WorkflowFailed)
      }
    case Event(CallFailed(call, failure), NoFailureMessage) =>
      persistStatus(call, ExecutionStatus.Failed)
      goto(WorkflowFailed) using FailureMessage(failure)
    case Event(Complete, NoFailureMessage) => goto(WorkflowSucceeded)
  }

  when(WorkflowFailed) {
    case Event(GetFailureMessage, msg: FailureMessage) =>
      sender() ! msg
      stay()
  }

  // FIXME the comment below appears misplaced?
  // We're supporting GetOutputs for all states, so there's nothing particular to do here
  when(WorkflowSucceeded)(FSM.NullFunction)

  whenUnhandled {
    case Event(GetOutputs, _) =>
      // FIXME There should be a more efficient way of getting final workflow outputs than getting every output
      // FIXME variable in the workflow including all intermediates between calls.
      val futureOutputs = for {
        symbols <- dataAccess.getOutputs(workflow.id)
      } yield (symbols map symbolStoreEntryToMapEntry).toMap
      futureOutputs pipeTo sender
      stay()
    case Event(e, _) =>
      log.debug(s"$tag received unhandled event $e while in state $stateName")
      stay()
  }

  // FSM will call *all* onTransition handlers which are defined for a particular state transition.
  // This handler will update workflow state for all transitions.
  onTransition {
    case fromState -> toState =>
      log.info(s"$tag transitioning from $fromState to $toState.")
      dataAccess.updateWorkflowState(workflow.id, toState)
  }

  private def persistStatus(call: Call, callStatus: CallStatus): Future[Unit] = {
    persistStatus(Iterable(call), callStatus)
  }

  private def persistStatus(calls: Traversable[Call], callStatus: CallStatus): Future[Unit] = {
    executionStore ++= calls.map { _ -> callStatus}.toMap
    val callFqns = calls.map(_.fullyQualifiedName)
    log.info(s"$tag persisting status of calls ${callFqns.mkString(", ")} to $callStatus.")
    dataAccess.setStatus(workflow.id, callFqns, callStatus)
  }

  private def awaitCallComplete(call: Call, outputs: CallOutputs): Try[Unit] = {
    val callFuture = handleCallCompleted(call, outputs)
    Await.ready(callFuture, DatabaseTimeout)
    callFuture.value.get
  }

  private def handleCallCompleted(call: Call, outputs: CallOutputs): Future[Unit] = {
    log.info(s"$tag handling completion of call '${call.fullyQualifiedName}'.")
    for {
      _ <- dataAccess.setOutputs(workflow.id, call, outputs)
      _ <- persistStatus(call, ExecutionStatus.Done)
    } yield ()
  }

  private def startActor(call: Call, locallyQualifiedInputs: Map[String, WdlValue]): Unit = {
    if (locallyQualifiedInputs.nonEmpty) {
      log.info(s"$tag inputs for call '${call.fullyQualifiedName}':")
      locallyQualifiedInputs foreach { input =>
        log.info(s"$tag     ${input._1} -> ${input._2}")
      }
    } else {
      log.info(s"$tag no inputs for call '${call.fullyQualifiedName}'")
    }

    val callActorProps = CallActor.props(call, locallyQualifiedInputs, backend, workflow)
    val callActor = context.actorOf(callActorProps)
    callActor ! CallActor.Start
    log.info(s"$tag created call actor for ${call.fullyQualifiedName}.")
  }

  /**
   * Start all calls which are currently in state `NotStarted` and whose prerequisites are all `Done`,
   * i.e. the calls which should now be eligible to run.  Update status to `Starting` in DataAccess and
   * start Call actors.  Return a Future of Iterable[Unit] to allow the caller to compose or block.
   */
  type CallAndInputs = (Call, Map[String, WdlValue])
  private def tryStartingRunnableCalls(): Try[Any] = {

    def isRunnable(call: Call) = call.prerequisiteCalls(workflow.namespace).forall(executionStore.get(_).get == ExecutionStatus.Done)
    /**
     * Start all calls which are currently in state `NotStarted` and whose prerequisites are all `Done`,
     * i.e. the calls which should now be eligible to run.
     */
    def findRunnableCalls: Iterable[Call] = {
      for {
        callEntry <- executionStore if callEntry._2 == ExecutionStatus.NotStarted
        call = callEntry._1 if isRunnable(call)
      } yield call
    }

    val runnableCalls = findRunnableCalls
    if (runnableCalls.isEmpty) {
      log.info(s"$tag no runnable calls to start.")
      Success(())
    } else {
      log.info(s"$tag starting calls: " + runnableCalls.map {_.fullyQualifiedName}.toSeq.sorted.mkString(", "))
      val futureCallsAndInputs = for {
        _ <- persistStatus(runnableCalls, ExecutionStatus.Starting)
        allInputs <- Future.sequence(runnableCalls map fetchLocallyQualifiedInputs)
      } yield runnableCalls zip allInputs
      Await.ready(futureCallsAndInputs, DatabaseTimeout)
      futureCallsAndInputs.value.get match {
        case Success(callsAndInputs) =>
          callsAndInputs foreach { case (call, inputs) => startActor(call, inputs) }
          Success(())
        case failure => failure
      }
    }
  }

  def fetchLocallyQualifiedInputs(call: Call): Future[Map[String, WdlValue]] = {
    /* TODO: This assumes that each Call has a parent that's a Workflow, but with scatter-gather that might not be the case */
    val parentWorkflow = call.parent.map {_.asInstanceOf[Workflow]} getOrElse {
      throw new WdlExpressionException("Expecting 'call' to have a 'workflow' parent.")
    }

    object CallInputWdlFunctions extends WdlFunctions {
      def getFunction(name: String): WdlFunction = {
        throw new WdlExpressionException("TODO: Some functions may be allowed in this context")
      }
    }

    def lookup(identifierString: String): WdlValue = {
      /* Start by assuming the identifierString is a Namespace reference, and if it is return that reference.
       * Otherwise, assume it's referencing a Call and return the value of that call as a WdlObject,
       * which represents the outputs of the Call (WdlObject only returned if outputs are present)
       */

      val namespaces = workflow.namespace.namespaces filter {_.importedAs.contains(identifierString)}
      namespaces.headOption.getOrElse {
        val matchedCall = parentWorkflow.calls.find {_.name == identifierString}.getOrElse {
          throw new WdlExpressionException(s"Expecting to find a call with name '$identifierString'")
        }
        val futureValue = fetchCallOutputEntries(matchedCall).map { entries =>
          val callOutputs = entries.map { entry =>
            val value = entry.wdlValue match {
              case Some(v) => v
              case _ => throw new WdlExpressionException(s"Could not evaluate call '${matchedCall.name}', because '${entry.key.name}' is undefined")
            }
            entry.key.name -> value
          }
          WdlObject(callOutputs.toMap)
        }
        Await.result(futureValue, DatabaseTimeout)
      }
    }

    fetchCallInputEntries(call).map { entries =>
      entries.map { entry =>
        val value = entry.wdlValue match {
          case Some(e: WdlExpression) => e.evaluate(lookup, CallInputWdlFunctions).get
          case Some(v) => v
          case _ => throw new WdlExpressionException("Unknown error")
        }
        entry.key.name -> value
      }.toMap
    }
  }

  private def fetchCallOutputEntries(call: Call) = dataAccess.getOutputs(workflow.id, call)
  private def fetchCallInputEntries(call: Call) = dataAccess.getInputs(workflow.id, call)

  /**
   * Load whatever execution statuses are stored for this workflow, regardless of whether this is a workflow being
   * restarted, or started for the first time.
   */
  private def createStore: Future[Map[Call, ExecutionStatus.Value]] = {
    dataAccess.getExecutionStatuses(workflow.id) map { statuses =>
      workflow.namespace.workflow.calls.map {call =>
        call -> statuses.get(call.fullyQualifiedName).get}.toMap
    }
  }

  private def buildSymbolStoreEntries(namespace: NamespaceWithWorkflow, inputs: HostInputs): Traversable[SymbolStoreEntry] = {
    val inputSymbols = inputs.map {case (name, value) => SymbolStoreEntry(name, value, input = true)}

    val callSymbols = for {
      call <- namespace.workflow.calls
      (k, v) <- call.inputMappings
    } yield SymbolStoreEntry(s"${call.fullyQualifiedName}.$k", v, input = true)

    inputSymbols.toSet ++ callSymbols.toSet
  }

  def createWorkflow(): Future[Unit] = {
    // TODO sanify WorkflowInfo/WorkflowDescriptor, the latter largely duplicates the former now.
    val workflowInfo = WorkflowInfo(workflow.id, workflow.wdlSource, workflow.wdlJson)
    // This only does the initialization for a newly created workflow.  For a restarted workflow we should be able
    // to assume the adjusted symbols already exist in the DB, but is it safe to assume the staged files are in place?
    val adjustedInputs = backend.initializeForWorkflow(workflow)
    dataAccess.createWorkflow(workflowInfo, buildSymbolStoreEntries(workflow.namespace, adjustedInputs), workflow.namespace.workflow.calls, backend)
  }

  private def isWorkflowDone: Boolean = executionStore.forall(_._2 == ExecutionStatus.Done)

  private def symbolStoreEntryToMapEntry(e: SymbolStoreEntry): (String, WdlValue) = {
    e.key.scope + "." + e.key.name -> e.wdlValue.get
  }
}
