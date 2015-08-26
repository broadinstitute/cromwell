package cromwell.engine.workflow

import akka.actor.{FSM, LoggingFSM, Props}
import akka.event.Logging
import akka.pattern.pipe
import cromwell.binding._
import cromwell.binding.types.WdlArrayType
import cromwell.binding.values.{WdlArray, WdlObject, WdlValue}
import cromwell.engine.ExecutionIndex._
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine._
import cromwell.engine.backend.Backend
import cromwell.engine.db.DataAccess.WorkflowInfo
import cromwell.engine.db.{CallStatus, DataAccess, ExecutionDatabaseKey}
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
  case class AbortComplete(call: OutputKey) extends WorkflowActorMessage
  case class CallStarted(call: OutputKey) extends WorkflowActorMessage
  case class CallCompleted(call: OutputKey, callOutputs: CallOutputs) extends WorkflowActorMessage
  case class CallFailed(call: OutputKey, failure: String) extends WorkflowActorMessage

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

  def isTerminal(es: ExecutionStatus): Boolean = TerminalStates contains es
  def isDone(entry: ExecutionStoreEntry) = entry._2 == ExecutionStatus.Done
  def isShard(key: CallKey): Boolean = key.index.isDefined

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
   * Try to generate output for a collector call, by collecting outputs for all of its shards.
   * It's fail-fast on shard output retrieval
   */
  private def generateCollectorOutput(collector: CollectorKey, shards: Iterable[CallKey]): Try[CallOutputs] = Try {
    val shardsOutputs = shards.toSeq sortBy { _.index.fromIndex } map { e =>
      fetchCallOutputEntries(e) map { _.value } get
    }
    collector.scope.task.outputs map { taskOutput =>
      val wdlValues = shardsOutputs.map(_.get(taskOutput.name).get)
      taskOutput.name -> new WdlArray(WdlArrayType(taskOutput.wdlType), wdlValues)
    } toMap
  }

  private def findShards(key: CollectorKey) = executionStore collect {
    case (k: CallKey, _) if k.scope == key.scope && isShard(k) => k
  }

  /**
   * Gathers and completes Collector Calls.
   */
  private def collectScatter(collector: CollectorKey) = {
    val shards: Iterable[CallKey] = findShards(collector)
    val collection = for {
      output <- generateCollectorOutput(collector, shards)
      callComplete <- awaitCallComplete(collector, output)
    } yield callComplete

    collection match {
      case Failure(e) =>
        self ! Event(CallFailed(collector, e.getMessage), NoFailureMessage)
      case Success(_) =>
        log.info(s"Collection complete for Scattered Call ${collector.scope.fullyQualifiedName}.")
    }
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
      executionsMarkdownTable().foreach {table => log.info(s"Initial executions:\n\n$table")}
      startRunnableCalls()
  }

  when(WorkflowRunning) {
    case Event(CallStarted(callKey), NoFailureMessage) =>
      persistStatus(callKey, ExecutionStatus.Running)
      stay()
    case Event(CallCompleted(callKey, outputs), NoFailureMessage) =>
      awaitCallComplete(callKey, outputs) match {
        case Success(_) =>
          if (isWorkflowDone) goto(WorkflowSucceeded) else startRunnableCalls()
        case Failure(e) =>
          log.error(e, e.getMessage)
          goto(WorkflowFailed)
      }
    case Event(CallFailed(callKey, failure), NoFailureMessage) =>
      persistStatus(callKey, ExecutionStatus.Failed)
      goto(WorkflowFailed) using FailureMessage(failure)
    case Event(Complete, NoFailureMessage) => goto(WorkflowSucceeded)
    case Event(AbortComplete(callKey), NoFailureMessage) =>
      // Something funky's going on if aborts are coming through while the workflow's still running. But don't second-guess
      // by transitioning the whole workflow - the message is either still in the queue or this command was maybe
      // cancelled by some external system.
      persistStatus(callKey, ExecutionStatus.Aborted)
      log.warning(s"Call ${callKey.scope.name} was aborted but the workflow should still be running.")
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
    case Event(AbortComplete(callKey), NoFailureMessage) =>
      persistStatus(callKey, ExecutionStatus.Aborted)
      if (isWorkflowAborted) goto(WorkflowAborted) using NoFailureMessage else stay()
    case Event(CallFailed(callKey, failure), NoFailureMessage) =>
      persistStatus(callKey, ExecutionStatus.Failed)
      if (isWorkflowAborted) goto(WorkflowAborted) using NoFailureMessage else stay()
    case Event(CallCompleted(callKey, outputs), NoFailureMessage) =>
      awaitCallComplete(callKey, outputs)
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

  private def persistStatus(key: OutputKey, callStatus: CallStatus): Future[Unit] = {
    persistStatus(Iterable(key), callStatus)
  }

  private def persistStatus(key: Traversable[OutputKey], callStatus: CallStatus): Future[Unit] = {
    executionStore ++= key map {
      _ -> callStatus
    }
    val callFqns = key.map(_.scope.fullyQualifiedName).mkString(", ")
    log.info(s"$tag persisting status of calls $callFqns to $callStatus.")
    dataAccess.setStatus(workflow.id, key map { k => ExecutionDatabaseKey(k.scope.fullyQualifiedName, k.index)}, callStatus)
  }

  private def awaitCallComplete(key: OutputKey, outputs: CallOutputs): Try[Unit] = {
    val callFuture = handleCallCompleted(key, outputs)
    Await.ready(callFuture, DatabaseTimeout)
    callFuture.value.get
  }

  private def handleCallCompleted(key: OutputKey, outputs: CallOutputs): Future[Unit] = {
    log.info(s"$tag handling completion of call '${key.scope.fullyQualifiedName}'.")
    for {
      _ <- dataAccess.setOutputs(workflow.id, key, outputs)
      _ <- persistStatus(key, ExecutionStatus.Done)
    } yield ()
  }

  private def startActor(callKey: CallKey, locallyQualifiedInputs: Map[String, WdlValue]): Unit = {
    if (locallyQualifiedInputs.nonEmpty) {
      log.info(s"$tag inputs for call '${callKey.scope.fullyQualifiedName}':")
      locallyQualifiedInputs foreach { input =>
        log.info(s"$tag     ${input._1} -> ${input._2}")
      }
    } else {
      log.info(s"$tag no inputs for call '${callKey.scope.fullyQualifiedName}'")
    }

    val callActorProps = CallActor.props(callKey, locallyQualifiedInputs, backend, workflow)
    val callActor = context.actorOf(callActorProps)
    callActor ! CallActor.Start
    log.info(s"$tag created call actor for ${callKey.scope.fullyQualifiedName}.")
  }

  private def tryStartingRunnableCalls(): Try[Any] = {

    def upstreamEntries(entry: ExecutionStoreKey, prerequisiteScope: Scope): Seq[ExecutionStoreEntry] = {
      prerequisiteScope.closestCommonAncestor(entry.scope) match {
        /**
         * If this entry refers to a Scope which has a common ancestor with prerequisiteScope
         * and that common ancestor is a Scatter block, then find the shard with the same index
         * as 'entry'.  In other words, if you're in the same scatter block as your pre-requisite
         * scope, then depend on the shard (with same index).
         *
         * NOTE: this algorithm was designed for ONE-LEVEL of scattering and probably does not
         * work as-is for nested scatter blocks
         */
        case Some(ancestor:Scatter) =>
          executionStore filter { case(k, v) => k.scope == prerequisiteScope && k.index == entry.index } toSeq

        /**
         * Otherwise, simply refer to the entry the collector entry.  This means that 'entry' depends
         * on the every shard of the pre-requisite scope to finish.
         */
        case _ =>
          executionStore filter { case(k, v) => k.scope == prerequisiteScope && k.index.isEmpty } toSeq
      }
    }

    def arePrerequisitesDone(key: ExecutionStoreKey): Boolean = {
      val upstream = key.scope.prerequisiteScopes.flatMap(s => upstreamEntries(key, s))
      upstream collect { case entry if entry._2 != ExecutionStatus.Done => entry } match {
        case s: Set[ExecutionStoreEntry] if s.nonEmpty => false
        case _ => true
      }
    }

    def isRunnable(entry: ExecutionStoreEntry) = {
      val (key, status) = entry
      status == ExecutionStatus.NotStarted && arePrerequisitesDone(key)
    }

    def findRunnableEntries: Traversable[ExecutionStoreEntry] = executionStore filter isRunnable

    val runnableEntries = findRunnableEntries

    val runnableCalls = runnableEntries collect { case(k: CallKey, v) => k.scope }
    if (runnableCalls.nonEmpty)
      log.info(s"$tag starting calls: " + runnableCalls.map(_.fullyQualifiedName).toSeq.sorted.mkString(", "))

    val entries = runnableEntries map {
      case (k: ScatterKey, v) =>
        /**
         * TODO (scatter) lookup function / WDL functions need to be real.  We're probably safe with
         * using NoFunctions, since a vast majority of the expression will be something like:
         * `x.y` -or- `x` -or- `[1,2,3]`, none of which are functions.  The lookup function can be crude
         * at first... simply interpreting variables as either references to Calls or Declarations.
         *
         * A more long-term approach would be to traverse the scope hierarchy to resolve a variable into
         * the closest definition in scope.
         */
        val rootWorkflow = k.scope.rootScope match {
          case w: Workflow => w
          case _ => throw new WdlExpressionException(s"Expected scatter '$k' to have a workflow root scope.")
        }

        def scatterLookupFunction(identifier: String): WdlValue = {
          resolveIdentifierOrElse(identifier, lookupCall(rootWorkflow), lookupDeclaration(rootWorkflow)) {
            throw new WdlExpressionException(s"Could not resolve identifier '$identifier' as a call or declaration.")
          }
        }

        val collection = k.scope.collection.evaluate(scatterLookupFunction, new NoFunctions )
        collection match {
          case Success(a: WdlArray) =>
            val newEntries = k.populate(a.value.size)
            persistStatus(newEntries.keys collect { case o: OutputKey => o }, ExecutionStatus.Starting)
            Success(newEntries)
          case Success(v: WdlValue) => Failure(new Throwable("Scatter collection must evaluate to an array"))
          case Failure(ex) => Failure(ex)
        }
      case (k: CollectorKey, v) =>
        // TODO (scatter) do collector stuff
        Failure(new UnsupportedOperationException("IMPLEMENT THIS"))
      case (k: CallKey, v) =>
        persistStatus(k, ExecutionStatus.Starting)
        val futureCallInputs = for {
          allInputs <- fetchLocallyQualifiedInputs(k.scope)
        } yield allInputs
        Await.ready(futureCallInputs, DatabaseTimeout)
        futureCallInputs.value.get match {
          case Success(callInputs) =>
            startActor(k, callInputs)
            Success(())
          case failure => failure
        }
      case (k, v) =>
        val message = s"$tag Unknown entry in execution store:\nKEY: $k\nVALUE:$v"
        log.error(message)
        Failure(new UnsupportedOperationException(message))
    }
    entries.find(_.isFailure) match {
      case Some(failure) => failure
      case _ => Success(())
    }
  }

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

  private def executionsMarkdownTable(): Option[String] = {
    val header = Seq("SCOPE", "INDEX", "TYPE", "STATUS")
    val rows = executionStore map { case (k, v) =>
      Seq(
        k.scope.fullyQualifiedNameWithIndexScopes.toString,
        k.index.getOrElse("").toString,
        k.getClass.getName,
        v.toString
      )
    }
    rows.toSeq match {
      case r:Seq[Seq[String]] if r.isEmpty => None
      case _ => Some(TerminalUtil.mdTable(rows.toSeq, header))
    }
  }

  private def lookupNamespace(name: String): Try[WdlNamespace] = {
    workflow.namespace.namespaces find { _.importedAs.contains(name) } match {
      case Some(x) => Success(x)
      case _ => Failure(new WdlExpressionException(s"Could not resolve $name as a namespace"))
    }
  }

  private def lookupCall(workflow: Workflow)(name: String): Try[WdlObject] = {
    workflow.calls find { _.name == name } match {
      case Some(matchedCall) =>
        /* FIXME: patch waiting for lookup function to be index aware */
        fetchCallOutputEntries(findCallKey(matchedCall, None).get)
      case None => Failure(new WdlExpressionException(s"Could not find a call with name '$name'"))
    }
  }

  private def lookupDeclaration(workflow: Workflow)(name: String): Try[WdlValue] = {
    workflow.declarations find { _.name == name } match {
      case Some(declaration) => fetchFullyQualifiedName(declaration.fullyQualifiedName)
      case None => Failure(new WdlExpressionException(s"Could not find a declaration with name '$name'"))
    }
  }

  private def resolveIdentifierOrElse(identifierString: String, resolvers: ((String) => Try[WdlValue]) *)(orElse: => Try[WdlValue]): WdlValue = {
    /* Try each of the resolver functions in order.  This uses a lazy Stream to only call a resolver function if a
     * preceding resolver function failed to resolve the identifier. */
    val attemptedResolutions = Stream(resolvers: _*) map { _(identifierString) } find { _.isSuccess }

    /* Return the first successful function's value or throw an exception. */
    attemptedResolutions.getOrElse {orElse}.get
  }

  private def findCallKey(call: Call, index: Option[Int]): Option[CallKey] = {
    executionStore find {
      case (k,v) => k.scope == call && k.index == index
    } collect { case (k: CallKey,v) => k }
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

    def lookup(identifier: String): WdlValue = {
      /* This algorithm defines three ways to lookup an identifier in order of their precedence:
       *
       *   1) Look for a WdlNamespace with matching name
       *   2) Look for a Call with a matching name (perhaps using a scope resolution algorithm)
       *   3) Look for a Declaration with a matching name (perhaps using a scope resolution algorithm)
       *
       * Each method is tried individually and the first to return a Success value takes precedence.
       */
      resolveIdentifierOrElse(identifier, lookupNamespace, lookupCall(parentWorkflow), lookupDeclaration(parentWorkflow)) {
        throw new WdlExpressionException(s"Could not resolve $identifier as a namespace, call, or declaration")
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

  private def fetchCallOutputEntries(callKey: CallKey): Try[WdlObject] = {
    val futureValue = dataAccess.getOutputs(workflow.id, ExecutionDatabaseKey(callKey.scope.fullyQualifiedName, callKey.index)).map {callOutputEntries =>
      val callOutputsAsMap = callOutputEntries.map { entry =>
        entry.key.name -> entry.wdlValue
      }.toMap
      callOutputsAsMap.find {case (k,v) => v.isEmpty} match {
        case Some(noneValue) => Failure(new WdlExpressionException(s"Could not evaluate call ${callKey.scope.name} because some if its inputs are not defined (i.e. ${noneValue._1}"))
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
        CallKey(call, None, None) -> statuses.get(ExecutionDatabaseKey(call.fullyQualifiedName, None)).get
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
    dataAccess.createWorkflow(workflowInfo, buildSymbolStoreEntries(workflow.namespace, adjustedInputs), workflow.namespace.workflow.children, backend)
  }

  private def isWorkflowDone: Boolean = executionStore forall isDone

  private def isWorkflowAborted: Boolean = executionStore forall { x => isTerminal(x._2) || x._2 == ExecutionStatus.NotStarted }
}
