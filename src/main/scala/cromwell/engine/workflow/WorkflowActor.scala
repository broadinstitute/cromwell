package cromwell.engine.workflow

import akka.actor.{FSM, LoggingFSM, Props}
import akka.event.Logging
import akka.pattern.pipe
import cromwell.binding._
import cromwell.binding.values.{WdlObject, WdlValue}
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine._
import cromwell.engine.backend.Backend
import cromwell.engine.db.DataAccess.WorkflowInfo
import cromwell.engine.db.{CallStatus, DataAccess}
import cromwell.engine.workflow.WorkflowActor._
import cromwell.util.TerminalUtil

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
  case object AbortWorkflow extends WorkflowActorMessage
  case class AbortComplete(call: Call) extends WorkflowActorMessage
  case class CallStarted(call: Call) extends WorkflowActorMessage
  case class CallCompleted(call: Call, callOutputs: CallOutputs) extends WorkflowActorMessage
  case class CallFailed(call: Call, failure: String) extends WorkflowActorMessage

  def props(descriptor: WorkflowDescriptor, backend: Backend, dataAccess: DataAccess): Props = {
    Props(WorkflowActor(descriptor, backend, dataAccess))
  }

  sealed trait WorkflowFailure
  case object NoFailureMessage extends WorkflowFailure
  case class FailureMessage(msg: String) extends WorkflowFailure with WorkflowActorMessage

  val DatabaseTimeout = 5 seconds

  type ExecutionStore = Map[ExecutionStoreKey, ExecutionStatus]
  type ExecutionStoreEntry = (ExecutionStoreKey, ExecutionStatus)

  // See also `WorkflowActor.createStore` that is similar, but uses the database.
  // Currently only used for testing. May be merged / replaced with createStore.
  def populate(workflow: Workflow): ExecutionStore = {
    (workflow.children map {
      case call: Call =>
        CallKey(call, None, None)
      case scatter: Scatter =>
        ScatterKey(scatter, None, None)
    } map {
      _ -> ExecutionStatus.NotStarted
    }).toMap
  }

  val TerminalStates = Vector(ExecutionStatus.Failed, ExecutionStatus.Done, ExecutionStatus.Aborted)
  def isExecutionStateFinished(es: ExecutionStatus): Boolean = {
    TerminalStates contains es
  }
}

case class WorkflowActor(workflow: WorkflowDescriptor,
                         backend: Backend,
                         dataAccess: DataAccess)
  extends LoggingFSM[WorkflowState, WorkflowFailure] with CromwellActor {
  private var executionStore: ExecutionStore = _
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
      symbolsMarkdownTable().foreach {table => log.info(s"Initial symbols:\n\n$table")}
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
    case Event(AbortComplete(call), NoFailureMessage) =>
      // Something funky's going on if aborts are coming through while the workflow's still running. But don't second-guess
      // by transitioning the whole workflow - the message is either still in the queue or this command was maybe
      // cancelled by some external system.
      persistStatus(call, ExecutionStatus.Aborted)
      log.warning(s"Call ${call.name} was aborted but the workflow should still be running.")
      stay()
  }

  when(WorkflowFailed) {
    case Event(GetFailureMessage, msg: FailureMessage) =>
      sender() ! msg
      stay()
  }

  // FIXME the comment below appears misplaced?
  // We're supporting GetOutputs for all states, so there's nothing particular to do here
  when(WorkflowSucceeded)(FSM.NullFunction)

  when(WorkflowAborting) {
    case Event(AbortComplete(call), NoFailureMessage) =>
      persistStatus(call, ExecutionStatus.Aborted)
      if (isWorkflowAborted) goto(WorkflowAborted) using NoFailureMessage else stay()
    case Event(CallFailed(call, failure), NoFailureMessage) =>
      persistStatus(call, ExecutionStatus.Failed)
      if (isWorkflowAborted) goto(WorkflowAborted) using NoFailureMessage else stay()
    case Event(CallCompleted(call, outputs), NoFailureMessage) =>
      awaitCallComplete(call, outputs)
      if (isWorkflowAborted) goto(WorkflowAborted) using NoFailureMessage else stay()
    case m =>
      log.error("Unexpected message in Aborting state: " + m.getClass.getSimpleName)
      if (isWorkflowAborted) goto(WorkflowAborted) using NoFailureMessage else stay()
  }

  when(WorkflowAborted)(FSM.NullFunction)

  whenUnhandled {
    case Event(GetOutputs, _) =>
      // NOTE: This is currently only used by tests
      // FIXME There should be a more efficient way of getting final workflow outputs than getting every output
      // FIXME variable in the workflow including all intermediates between calls.
      val futureOutputs = dataAccess.getOutputs(workflow.id) map SymbolStoreEntry.toWorkflowOutputs
      futureOutputs pipeTo sender
      stay()
    case Event(AbortWorkflow, _) =>
      context.children foreach { _ ! CallActor.AbortCall }
      goto(WorkflowAborting) using NoFailureMessage
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
    executionStore ++= calls map {
      CallKey(_, None, None) -> callStatus
    }
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
    val callActor = context.actorOf(callActorProps, s"CallActor-${workflow.id}-${call.name}")
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

    def isRunnable(scope: Scope) = {
      scope match {
        case call: Call =>
          // FIXME: Sticking w/ prerequisiteCalls for now until our system is more scope aware
          call.prerequisiteCalls forall {
            case previousCall =>
              executionStore.find(_._1.scope == previousCall).get._2 == ExecutionStatus.Done
          }
      }
    }
    /**
     * Start all calls which are currently in state `NotStarted` and whose prerequisites are all `Done`,
     * i.e. the calls which should now be eligible to run.
     */
    def findRunnableCalls: Iterable[Call] = {
      for {
      // TODO: Need other scopes
        callEntry <- executionStore if callEntry._1.scope.isInstanceOf[Call] && callEntry._2 == ExecutionStatus.NotStarted
        call = callEntry._1.scope.asInstanceOf[Call] if isRunnable(call)
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
  /* Return a Markdown table of all entries in the database */
  private def symbolsMarkdownTable(): Option[String] = {
    val header = Seq("SCOPE", "NAME", "I/O", "TYPE", "VALUE")
    val maxColChars = 100
    val rows = fetchAllEntries.map { entry =>
      val valueString = entry.wdlValue match {
        case Some(value) => s"(${value.wdlType.toWdlString}) " + value.valueString
        case _ => ""
      }
      Seq(
        entry.key.scope,
        entry.key.name,
        if (entry.key.input) "INPUT" else "OUTPUT",
        entry.wdlType.toWdlString,
        if (valueString.length > maxColChars) valueString.substring(0, maxColChars) else valueString
      )
    }.toSeq

    rows match {
      case r: Seq[Seq[String]] if r.isEmpty => None
      case _ => Some(TerminalUtil.mdTable(rows, header))
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
      /* This algorithm defines three ways to lookup an identifier in order of their precedence:
       *
       *   1) Look for a WdlNamespace with matching name
       *   2) Look for a Call with a matching name (perhaps using a scope resolution algorithm)
       *   3) Look for a Declaration with a matching name (perhaps using a scope resolution algorithm)
       *
       * Each method is tried individually and the first to return a Success value takes precedence.
       */

      def lookupNamespace(name: String): Try[WdlNamespace] = {
        workflow.namespace.namespaces.find {_.importedAs.contains(name)} match {
          case Some(x) => Success(x)
          case _ => Failure(new WdlExpressionException(s"Could not resolve $identifierString as a namespace"))
        }
      }

      def lookupCall(name: String): Try[WdlObject] = {
        parentWorkflow.calls.find {_.name == identifierString} match {
          case Some(matchedCall) => fetchCallOutputEntries(matchedCall)
          case None => Failure(new WdlExpressionException(s"Could not find a call with name '$identifierString'"))
        }
      }

      def lookupDeclaration(name: String): Try[WdlValue] = {
        parentWorkflow.declarations.find {_.name == identifierString} match {
          case Some(declaration) => fetchFullyQualifiedName(declaration.fullyQualifiedName)
          case None => Failure(new WdlExpressionException(s"Could not find a declaration with name '$identifierString'"))
        }
      }

      def notFound = throw new WdlExpressionException(s"Could not resolve $identifierString as a namespace, call, or declaration")

      /* Try each of the three methods in sequence. */
      val attemptedResolutions = Stream(lookupNamespace _, lookupCall _, lookupDeclaration _).map {_(identifierString)}.find {_.isSuccess}

      /* Return the first successful function's value or throw exception */
      attemptedResolutions.getOrElse {notFound}.get
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

  private def fetchFullyQualifiedName(fqn: FullyQualifiedName): Try[WdlValue] = {
    val futureValue = dataAccess.getFullyQualifiedName(workflow.id, fqn).map {
      case t: Traversable[SymbolStoreEntry] if t.isEmpty =>
        Failure(new WdlExpressionException(s"Could not find a declaration with fully-qualified name '$fqn'"))
      case t: Traversable[SymbolStoreEntry] if t.size > 1 =>
        Failure(new WdlExpressionException(s"Expected only one declaration for fully-qualified name '$fqn', got ${t.size}"))
      case t: Traversable[SymbolStoreEntry] => t.head.wdlValue match {
        case Some(value) => Success(value)
        case None => Failure(new WdlExpressionException(s"No value defined for fully-qualified name $fqn"))
      }
    }
    Await.result(futureValue, DatabaseTimeout)
  }

  private def fetchAllEntries: Traversable[SymbolStoreEntry] = {
    val futureValue = dataAccess.getAll(workflow.id)
    Await.result(futureValue, DatabaseTimeout)
  }

  private def fetchCallOutputEntries(call: Call): Try[WdlObject] = {
    val futureValue = dataAccess.getOutputs(workflow.id, call.fullyQualifiedName).map {callOutputEntries =>
      val callOutputsAsMap = callOutputEntries.map { entry =>
        entry.key.name -> entry.wdlValue
      }.toMap
      callOutputsAsMap.find {case (k,v) => v.isEmpty} match {
        case Some(noneValue) => Failure(new WdlExpressionException(s"Could not evaluate call ${call.name} because some if its inputs are not defined (i.e. ${noneValue._1}"))
        case None => Success(WdlObject(callOutputsAsMap.map {case (k,v) => k -> v.get}))
      }
    }
    Await.result(futureValue, DatabaseTimeout)
  }

  private def fetchCallInputEntries(call: Call) = dataAccess.getInputs(workflow.id, call)

  /**
   * Load whatever execution statuses are stored for this workflow, regardless of whether this is a workflow being
   * restarted, or started for the first time.
   */
  private def createStore: Future[ExecutionStore] = {
    dataAccess.getExecutionStatuses(workflow.id) map { statuses =>
      // TODO: This is NOT right. For example: statuses from DB are not just for a FQN, need the scatter index.
      // See also `WorkflowActor.populate` that is similar, but uses indexes.
      workflow.namespace.workflow.calls.map(call =>
        CallKey(call, None, None) -> statuses.get(call.fullyQualifiedName).get
      ).toMap
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

  private def isWorkflowAborted: Boolean = executionStore forall { x => isExecutionStateFinished(x._2) || x._2 == ExecutionStatus.NotStarted }
}
