package cromwell.engine.workflow

import akka.actor.{FSM, LoggingFSM, Props}
import akka.event.Logging
import cromwell.binding._
import cromwell.binding.expression.NoFunctions
import cromwell.binding.types.WdlArrayType
import cromwell.binding.values.{WdlArray, WdlCallOutputsObject, WdlValue}
import cromwell.engine.ExecutionIndex._
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine._
import cromwell.engine.backend.Backend
import cromwell.engine.db.DataAccess._
import cromwell.engine.db.{CallStatus, ExecutionDatabaseKey}
import cromwell.engine.workflow.WorkflowActor._
import cromwell.util.TerminalUtil

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration._
import scala.concurrent.{Await, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object WorkflowActor {
  sealed trait WorkflowActorMessage
  case object Complete extends WorkflowActorMessage
  case object GetFailureMessage extends WorkflowActorMessage
  case object AbortWorkflow extends WorkflowActorMessage
  case class AbortComplete(call: OutputKey) extends WorkflowActorMessage
  case class CallStarted(call: OutputKey) extends WorkflowActorMessage
  case class CallCompleted(call: OutputKey, callOutputs: CallOutputs, returnCode: Int) extends WorkflowActorMessage
  case class CallFailed(call: OutputKey, returnCode: Option[Int], failure: String) extends WorkflowActorMessage
  case object Terminate extends WorkflowActorMessage
  case class ExecutionStoreCreated(startMode: StartMode) extends WorkflowActorMessage
  case class AsyncFailure(t: Throwable) extends WorkflowActorMessage
  final case class PerformTransition(toState: WorkflowState) extends WorkflowActorMessage

  sealed trait StartMode {
    def runInitialization(actor: WorkflowActor): Future[Unit]
  }

  case object Start extends WorkflowActorMessage with StartMode {
    override def runInitialization(actor: WorkflowActor): Future[Unit] = {
      // This only does the initialization for a newly created workflow.  For a restarted workflow we should be able
      // to assume the adjusted symbols already exist in the DB, but is it safe to assume the staged files are in place?
      actor.initializeWorkflow match {
        case Success(inputs) => actor.createWorkflow(inputs)
        case Failure(ex) => Future.failed(ex)
      }
    }
  }
  
  case object Restart extends WorkflowActorMessage with StartMode {
    override def runInitialization(actor: WorkflowActor): Future[Unit] = Future.successful(())
  }

  def props(descriptor: WorkflowDescriptor, backend: Backend): Props = {
    Props(WorkflowActor(descriptor, backend))
  }

  sealed trait WorkflowFailure
  case object NoFailureMessage extends WorkflowFailure
  case class FailureMessage(msg: String) extends WorkflowFailure with WorkflowActorMessage

  val AkkaTimeout = 5 seconds

  type ExecutionStore = Map[ExecutionStoreKey, ExecutionStatus]
  type ExecutionStoreEntry = (ExecutionStoreKey, ExecutionStatus)

  val TerminalStates = Vector(ExecutionStatus.Failed, ExecutionStatus.Done, ExecutionStatus.Aborted)

  def isExecutionStateFinished(es: ExecutionStatus): Boolean = TerminalStates contains es

  val DontRepeatTimer = false

  def isTerminal(status: ExecutionStatus): Boolean = TerminalStates contains status
  def isDone(entry: ExecutionStoreEntry): Boolean = entry._2 == ExecutionStatus.Done
  def isShard(key: CallKey): Boolean = key.index.isDefined

  private val MarkdownMaxColumnChars = 100
}

case class WorkflowActor(workflow: WorkflowDescriptor,
                         backend: Backend)
  extends LoggingFSM[WorkflowState, WorkflowFailure] with CromwellActor {

  def createWorkflow(inputs: HostInputs): Future[Unit] = {
    globalDataAccess.createWorkflow(
      workflow, buildSymbolStoreEntries(workflow.namespace, inputs), workflow.namespace.workflow.children, backend)
  }

  private var executionStore: ExecutionStore = _
  val tag: String = s"WorkflowActor [UUID(${workflow.shortId})]"
  override val log = Logging(context.system, classOf[WorkflowActor])

  startWith(WorkflowSubmitted, NoFailureMessage)

 /**
   * Try to generate output for a collector call, by collecting outputs for all of its shards.
   * It's fail-fast on shard output retrieval
   */
  private def generateCollectorOutput(collector: CollectorKey, shards: Iterable[CallKey]): Try[CallOutputs] = Try {
    val shardsOutputs = shards.toSeq sortBy { _.index.fromIndex } map { e =>
      fetchCallOutputEntries(e) map { _.outputs } getOrElse(throw new RuntimeException(s"Could not retrieve output for shard ${e.scope} #${e.index}"))
    }
    collector.scope.task.outputs map { taskOutput =>
      val wdlValues = shardsOutputs.map(s => s.getOrElse(taskOutput.name, throw new RuntimeException(s"Could not retrieve output ${taskOutput.name}")))
      taskOutput.name -> new WdlArray(WdlArrayType(taskOutput.wdlType), wdlValues)
    } toMap
  }

  private def findShardEntries(key: CollectorKey): Iterable[ExecutionStoreEntry] = executionStore collect {
    case (k: CallKey, v) if k.scope == key.scope && isShard(k) => (k, v)
  }

  /**
   * Attempt to start all runnable calls and return the next FSM state.  If successful this will be
   * `WorkflowRunning`, otherwise `WorkflowFailed`.
   */
  private def startRunnableCalls(): State = {
    tryStartingRunnableCalls() match {
      case Success(entries) => if (entries.nonEmpty) startRunnableCalls() else goto(WorkflowRunning)
      case Failure(e) =>
        log.error(e, e.getMessage)
        goto(WorkflowFailed)
    }
  }

  private def initializeExecutionStore(initializationCode: Future[Unit], startMode: StartMode): Unit = {
    val futureStore = for {
      _ <- initializationCode
      store <- createStore
    } yield store

    futureStore onComplete {
      case Success(store) =>
        executionStore = store
        self ! ExecutionStoreCreated(startMode)
      case Failure(t) =>
        self ! AsyncFailure(t)
    }
  }

  private def initializeWorkflow: Try[HostInputs] = backend.initializeForWorkflow(workflow)

  /**
   * Dump symbol and execution tables, start runnable calls, and message self to transition to the appropriate
   * next state.
   */
  private def handleRestart(): Unit = {
    val workflowState = for {
      symbols <- symbolsMarkdownTable
      _ = symbols foreach { table => log.info(s"Initial symbols:\n\n$table") }
      executions <- executionsMarkdownTable
      _ = executions foreach { table => log.info(s"Initial executions:\n\n$table") }
      state = startRunnableCalls()
    } yield state.stateName

    workflowState onComplete {
      case Success(state) => self ! PerformTransition(state)
      case Failure(t) => self ! AsyncFailure(t)
    }
  }

  when(WorkflowSubmitted) {
    case Event(startMode: StartMode, NoFailureMessage) =>
      log.info(s"$tag $startMode message received")
      initializeExecutionStore(startMode.runInitialization(this), startMode)
      stay()
    case Event(ExecutionStoreCreated(startMode), NoFailureMessage) =>
      if (startMode == Restart) {
        handleRestart()
        stay()
      } else {
        startRunnableCalls()
      }
    case Event(PerformTransition(toState), NoFailureMessage) =>
      goto(toState)
    case Event(AsyncFailure(t), NoFailureMessage) =>
      log.error(t.getMessage, t)
      goto(WorkflowFailed)
  }

  when(WorkflowRunning) {
    case Event(CallStarted(callKey), NoFailureMessage) =>
      persistStatus(callKey, ExecutionStatus.Running)
      stay()
    case Event(CallCompleted(callKey, outputs, returnCode), NoFailureMessage) =>
      awaitCallComplete(callKey, outputs, returnCode) match {
        case Success(_) =>
          if (isWorkflowDone) goto(WorkflowSucceeded) else startRunnableCalls()
        case Failure(e) =>
          log.error(e, e.getMessage)
          goto(WorkflowFailed)
      }
    case Event(CallFailed(callKey, returnCode, failure), NoFailureMessage) =>
      persistStatus(callKey, ExecutionStatus.Failed, returnCode)
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

  when(WorkflowSucceeded) {
    case Event(Terminate, _) =>
      log.debug(s"$tag WorkflowActor is done, shutting down.")
      context.stop(self)
      stay()
  }

  when(WorkflowAborting) {
    case Event(AbortComplete(callKey), NoFailureMessage) =>
      persistStatus(callKey, ExecutionStatus.Aborted, None)
      if (isWorkflowAborted) goto(WorkflowAborted) using NoFailureMessage else stay()
    case Event(CallFailed(callKey, returnCode, failure), NoFailureMessage) =>
      persistStatus(callKey, ExecutionStatus.Failed, returnCode)
      if (isWorkflowAborted) goto(WorkflowAborted) using NoFailureMessage else stay()
    case Event(CallCompleted(callKey, outputs, returnCode), NoFailureMessage) =>
      awaitCallComplete(callKey, outputs, returnCode)
      if (isWorkflowAborted) goto(WorkflowAborted) using NoFailureMessage else stay()
    case m =>
      log.error("Unexpected message in Aborting state: " + m.getClass.getSimpleName)
      if (isWorkflowAborted) goto(WorkflowAborted) using NoFailureMessage else stay()
  }

  when(WorkflowAborted)(FSM.NullFunction)

  whenUnhandled {
    case Event(AbortWorkflow, _) =>
      context.children foreach { _ ! CallActor.AbortCall }
      goto(WorkflowAborting) using NoFailureMessage
    case Event(e, _) =>
      log.debug(s"$tag received unhandled event $e while in state $stateName")
      stay()
  }

  /*
    FSM will call *all* onTransition handlers which are defined for a particular state transition.
    This handler will update workflow state for all transitions.
  */
  onTransition {
    case fromState -> toState =>
      def handleTerminalWorkflow: Future[Unit] = {
        for {
          _ <- backend.cleanUpForWorkflow(workflow)
          _ <- globalDataAccess.updateWorkflowOptions(workflow.id, workflow.workflowOptions.clearEncryptedValues)
          //  Send a message to self to trigger an actor shutdown. Run on a short timer to help enable some
          //  unit test instrumentation
          _ = setTimer(s"WorkflowActor termination message: $tag", Terminate, AkkaTimeout, DontRepeatTimer)
        } yield ()
      }

      val transitionFuture = for {
        // Write the new workflow state before logging the change, tests assume the change is in effect when
        // the message is logged.
        _ <- globalDataAccess.updateWorkflowState(workflow.id, toState)
        _ = log.info(s"$tag transitioning from $fromState to $toState.")
        _ <- if (toState.isTerminal) handleTerminalWorkflow else Future.successful({})
      } yield ()

      transitionFuture recover {
        case e: Exception => log.error(e, s"Failed to transition workflow status from $fromState to $toState for $tag")
      }
  }

  private def persistStatus(key: ExecutionStoreKey, status: ExecutionStatus,
                            returnCode: Option[Int] = None): Unit = {
    Await.result(persistStatuses(Iterable(key), status, returnCode), Duration.Inf)
  }

  private def persistStatuses(key: Traversable[ExecutionStoreKey], executionStatus: ExecutionStatus,
                              returnCode: Option[Int] = None): Future[Unit] = {

    val databaseKeys = key map { k =>
      log.info(s"$tag persisting status of ${k.tag} to $executionStatus.")
      ExecutionDatabaseKey(k.scope.fullyQualifiedName, k.index)
    }

    for {
      // Write the status to the database before updating the store, the store is what is examined to
      // determine workflow doneness and if that persisted workflow representation is not consistent,
      // tests may see unexpected values.
      _ <- globalDataAccess.setStatus(workflow.id, databaseKeys, CallStatus(executionStatus, returnCode))
      _ = executionStore ++= key map { _ -> executionStatus }
    } yield ()
  }

  private def awaitCallComplete(key: OutputKey, outputs: CallOutputs, returnCode: Int): Try[Unit] = {
    val callFuture = handleCallCompleted(key, outputs, returnCode)
    Await.ready(callFuture, AkkaTimeout)
    callFuture.value.get
  }

  private def handleCallCompleted(key: OutputKey, outputs: CallOutputs, returnCode: Int): Future[Unit] = {
    log.info(s"$tag handling completion of call '${key.tag}'.")
    for {
      _ <- globalDataAccess.setOutputs(workflow.id, key, outputs)
      _ = persistStatus(key, ExecutionStatus.Done, Option(returnCode))
    } yield()
  }

  private def startActor(callKey: CallKey, locallyQualifiedInputs: Map[String, WdlValue]): Unit = {
    if (locallyQualifiedInputs.nonEmpty) {
      val inputs = locallyQualifiedInputs map { case(lqn, value) => s"  $lqn -> $value" } mkString "\n"
      log.info(s"$tag inputs for call '${callKey.tag}':\n$inputs")
    } else {
      log.info(s"$tag no inputs for call '${callKey.tag}'")
    }

    val callActorProps = CallActor.props(callKey, locallyQualifiedInputs, backend, workflow)
    val callActor = context.actorOf(callActorProps)
    callActor ! CallActor.Start
    log.info(s"$tag created call actor for ${callKey.tag}.")
  }

  private def tryStartingRunnableCalls(): Try[Traversable[ExecutionStoreKey]] = {

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
        case Some(ancestor: Scatter) =>
          executionStore filter { case(k, _) => k.scope == prerequisiteScope && k.index == entry.index } toSeq

        /**
         * Otherwise, simply refer to the entry the collector entry.  This means that 'entry' depends
         * on every shard of the pre-requisite scope to finish.
         */
        case _ =>
          executionStore filter { case(k, _) => k.scope == prerequisiteScope && k.index.isEmpty } toSeq
      }
    }

    def arePrerequisitesDone(key: ExecutionStoreKey): Boolean = {
      val upstream = key.scope.prerequisiteScopes.map(s => upstreamEntries(key, s))
      val downstream = key match {
        case collector: CollectorKey => findShardEntries(collector)
        case _ => Nil
      }
      val dependencies = upstream.flatten ++ downstream
      val dependenciesResolved = dependencies.isEmpty || dependencies.forall(isDone)

      /**
       * We need to make sure that all prerequisiteScopes have been resolved to some entry before going forward.
       * If a scope cannot be resolved it may be because it is in a scatter that has not been populated yet,
       * therefore there is no entry in the executionStore for this scope.
       * If that's the case this prerequisiteScope has not been run yet, hence the (upstream forall {_.nonEmpty})
       */
      (upstream forall { _.nonEmpty }) && dependenciesResolved
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

    val entries: Traversable[Try[Iterable[ExecutionStoreKey]]] = runnableEntries map {
      case (k: ScatterKey, v) => processRunnableScatter(k, v)
      case (k: CollectorKey, v) => processRunnableCollector(k)
      case (k: CallKey, v) => processRunnableCall(k, v)
      case (k, v) =>
        val message = s"$tag Unknown entry in execution store:\nKEY: $k\nVALUE:$v"
        log.error(message)
        Failure(new UnsupportedOperationException(message))
    }

    entries.find(_.isFailure) match {
      case Some(failure) => failure
      case _ => Success(entries.flatMap(_.get))
    }
  }

  private def lookupNamespace(name: String): Try[WdlNamespace] = {
    workflow.namespace.namespaces find { _.importedAs.contains(name) } match {
      case Some(x) => Success(x)
      case _ => Failure(new WdlExpressionException(s"Could not resolve $name as a namespace"))
    }
  }

  private def lookupCall(key: ExecutionStoreKey, workflow: Workflow)(name: String): Try[WdlCallOutputsObject] = {
    workflow.calls find { _.name == name } match {
      case Some(matchedCall) =>
        /**
         * After matching the Call, this determines if the `key` depends on a single shard
         * of a scatter'd job or if it depends on the whole thing.  Right now, the heuristic
         * is "If we're both in a scatter block together, then I depend on a shard.  If not,
         * I depend on the collected value"
         *
         * TODO: nested-scatter - this will likely not be sufficient for nested scatters
         */
        val index: ExecutionIndex = matchedCall.closestCommonAncestor(key.scope) flatMap {
          case s: Scatter => key.index
          case _ => None
        }
        fetchCallOutputEntries(findCallKey(matchedCall, index) getOrElse {
          throw new WdlExpressionException(s"Could not find a callKey for name '${matchedCall.name}'")
        })
      case None => Failure(new WdlExpressionException(s"Could not find a call with name '$name'"))
    }
  }

  private def lookupDeclaration(workflow: Workflow)(name: String): Try[WdlValue] = {
    workflow.declarations find { _.name == name } match {
      case Some(declaration) => fetchFullyQualifiedName(declaration.fullyQualifiedName)
      case None => Failure(new WdlExpressionException(s"Could not find a declaration with name '$name'"))
    }
  }

  private def lookupScatterVariable(callKey: CallKey, workflow: Workflow)(name: String): Try[WdlValue] = {
    val scatterBlock = callKey.scope.ancestry collect { case s: Scatter => s } find { _.item == name }
    val scatterCollection = scatterBlock map { s =>
      s.collection.evaluate(scatterCollectionLookupFunction(workflow, callKey), new NoFunctions) match {
        case Success(v: WdlArray) if callKey.index.isDefined =>
          if (v.value.isDefinedAt(callKey.index.get))
            Success(v.value(callKey.index.get))
          else
            Failure(new WdlExpressionException(s"Index ${callKey.index.get} out of bounds for $name array."))
        case Success(v: WdlArray) => Failure(new WdlExpressionException(s"$name evaluated to an Array but $callKey has no index"))
        case _ => Failure(new WdlExpressionException(s"$name did not evaluate to a WdlArray"))
      }
    }
    scatterCollection.getOrElse(
      Failure(new WdlExpressionException(s"$name is does not reference a scattered variable"))
    )
  }

  private def resolveIdentifierOrElse(identifierString: String, resolvers: ((String) => Try[WdlValue]) *)(orElse: => Try[WdlValue]): WdlValue = {
    /* Try each of the resolver functions in order.  This uses a lazy Stream to only call a resolver function if a
     * preceding resolver function failed to resolve the identifier. */
    val attemptedResolutions = Stream(resolvers: _*) map { _(identifierString) } find { _.isSuccess }

    /* Return the first successful function's value or throw an exception. */
    attemptedResolutions.getOrElse(orElse).get
  }

  private def findCallKey(call: Call, index: ExecutionIndex): Option[OutputKey] = {
    executionStore.collect({
      case (k: OutputKey, _) if k.scope == call && k.index == index => k
    }).headOption
  }

  def fetchLocallyQualifiedInputs(callKey: CallKey): Future[Map[String, WdlValue]] = {
    val parentWorkflow = callKey.scope.ancestry.lastOption map { _.asInstanceOf[Workflow] } getOrElse {
      throw new WdlExpressionException("Expecting 'call' to have a 'workflow' parent.")
    }

    def lookup(identifier: String): WdlValue = {
      /* This algorithm defines three ways to lookup an identifier in order of their precedence:
       *
       *   1) Traverse up the scope hierarchy and see if the variable reference any scatter item
       *   2) Look for a WdlNamespace with matching name
       *   3) Look for a Call with a matching name (perhaps using a scope resolution algorithm)
       *   4) Look for a Declaration with a matching name (perhaps using a scope resolution algorithm)
       *
       *  Each method is tried individually and the first to return a Success value takes precedence.
       */

      resolveIdentifierOrElse(identifier, lookupScatterVariable(callKey, parentWorkflow), lookupNamespace, lookupCall(callKey, parentWorkflow), lookupDeclaration(parentWorkflow)) {
        throw new WdlExpressionException(s"Could not resolve $identifier as a scatter variable, namespace, call, or declaration")
      }
    }

    fetchCallInputEntries(callKey.scope).map { entries =>
      entries.map { entry =>
        // .get are used below because the exception will be captured by the Future
        val taskInput = findTaskInput(entry.scope, entry.key.name).get

        val value = entry.wdlValue match {
          case Some(e: WdlExpression) => e.evaluate(lookup, new NoFunctions)
          case Some(v) => Success(v)
          case _ => Failure(new WdlExpressionException("Unknown error"))
        }
        val coercedValue = value.flatMap(x => taskInput.wdlType.coerceRawValue(x))
        entry.key.name -> coercedValue.get
      }.toMap
    }
  }

  private def findTaskInput(callFqn: String, inputName: String): Try[TaskInput] = {
    val exception = new WdlExpressionException(s"Could not find task input '$inputName' for call '$callFqn'")
    workflow.namespace.resolve(callFqn) match {
      case Some(c:Call) =>
        c.task.inputs.find(_.name == inputName) match {
          case Some(i) => Success(i)
          case None => Failure(exception)
        }
      case _ => Failure(exception)
    }
  }

  private def fetchFullyQualifiedName(fqn: FullyQualifiedName): Try[WdlValue] = {
    val futureValue = globalDataAccess.getFullyQualifiedName(workflow.id, fqn).map {
      case t: Traversable[SymbolStoreEntry] if t.isEmpty =>
        Failure(new WdlExpressionException(s"Could not find a declaration with fully-qualified name '$fqn'"))
      case t: Traversable[SymbolStoreEntry] if t.size > 1 =>
        Failure(new WdlExpressionException(s"Expected only one declaration for fully-qualified name '$fqn', got ${t.size}"))
      case t: Traversable[SymbolStoreEntry] => t.head.wdlValue match {
        case Some(value) => Success(value)
        case None => Failure(new WdlExpressionException(s"No value defined for fully-qualified name $fqn"))
      }
    }
    Await.result(futureValue, AkkaTimeout)
  }

  private def fetchCallOutputEntries(outputKey: OutputKey): Try[WdlCallOutputsObject] = {
    val futureValue = globalDataAccess.getOutputs(workflow.id, ExecutionDatabaseKey(outputKey.scope.fullyQualifiedName, outputKey.index)).map {callOutputEntries =>
      val callOutputsAsMap = callOutputEntries.map(entry => entry.key.name -> entry.wdlValue).toMap
      callOutputsAsMap find { case (k, v) => v.isEmpty } match {
        case Some(noneValue) => Failure(new WdlExpressionException(s"Could not evaluate call ${outputKey.scope.name} because some of its inputs are not defined (i.e. ${noneValue._1}"))
        // TODO: .asInstanceOf[Call]?
        case None => Success(WdlCallOutputsObject(outputKey.scope.asInstanceOf[Call], callOutputsAsMap.map {
          case (k, v) => k -> v.getOrElse {
            throw new WdlExpressionException(s"Could not retrieve output $k value for call ${outputKey.scope.name}")
          }
        }))
      }
    }
    Await.result(futureValue, AkkaTimeout)
  }

  private def fetchCallInputEntries(call: Call) = globalDataAccess.getInputs(workflow.id, call)

  /**
   * Load whatever execution statuses are stored for this workflow, regardless of whether this is a workflow being
   * restarted, or started for the first time.
   */
  private def createStore: Future[ExecutionStore] = {
    def isInScatterBlock(c: Call) = c.ancestry.exists(_.isInstanceOf[Scatter])
    globalDataAccess.getExecutionStatuses(workflow.id) map { statuses =>
      statuses map { case (k, v) =>
        val key: ExecutionStoreKey = (workflow.namespace.resolve(k.fqn), k.index) match {
          case (Some(c: Call), Some(i)) => CallKey(c, Option(i))
          case (Some(c: Call), None) if isInScatterBlock(c) => CollectorKey(c)
          case (Some(c: Call), None) => CallKey(c, None)
          case (Some(s: Scatter), None) => ScatterKey(s, None)
          case _ => throw new UnsupportedOperationException(s"Execution entry invalid: $k -> $v")
        }
        key -> v.executionStatus
      }
    }
  }

  private def buildSymbolStoreEntries(namespace: NamespaceWithWorkflow, inputs: HostInputs): Traversable[SymbolStoreEntry] = {
    val inputSymbols = inputs map { case (name, value) => SymbolStoreEntry(name, value, input = true) }

    val callSymbols = for {
      call <- namespace.workflow.calls
      (k, v) <- call.inputMappings
    } yield SymbolStoreEntry(s"${call.fullyQualifiedName}.$k", v, input = true)

    inputSymbols.toSet ++ callSymbols.toSet
  }

  /**
   * This is the lookup function used to evaluate scatter collection expressions.
   *
   * For example, scatter(x in foo.bar) would evaluate the collection "foo.bar"
   * and call this lookup function on "foo".
   *
   * This implementation takes a few shortcuts and tries to find a Call or
   * Declaration with the given name in the workflow.
   *
   * A more long-term approach would be to traverse the scope hierarchy to resolve a variable into
   * the closest definition in scope.
   */
  private def scatterCollectionLookupFunction(workflow: Workflow, key: ExecutionStoreKey)(identifier: String): WdlValue = {
    resolveIdentifierOrElse(identifier, lookupCall(key, workflow), lookupDeclaration(workflow)) {
      throw new WdlExpressionException(s"Could not resolve identifier '$identifier' as a call or declaration.")
    }
  }

  private def isWorkflowDone: Boolean = executionStore forall isDone

  private def isWorkflowAborted: Boolean = executionStore.values forall { state => isTerminal(state) || state == ExecutionStatus.NotStarted }

  private def processRunnableScatter(scatterKey: ScatterKey, status: ExecutionStatus): Try[Iterable[ExecutionStoreKey]] = {
    val rootWorkflow = scatterKey.scope.rootScope match {
      case w: Workflow => w
      case _ => throw new WdlExpressionException(s"Expected scatter '$scatterKey' to have a workflow root scope.")
    }

    val collection = scatterKey.scope.collection.evaluate(scatterCollectionLookupFunction(rootWorkflow, scatterKey), new NoFunctions)
    collection match {
      case Success(a: WdlArray) => Try {
        val newEntries = scatterKey.populate(a.value.size)
        persistStatus(scatterKey, ExecutionStatus.Starting, None)
        val createScatter = for {
          _ <- globalDataAccess.insertCalls(workflow.id, newEntries.keys, backend)
          _ <- persistStatuses(newEntries.keys, ExecutionStatus.NotStarted, None)
          _ = persistStatus(scatterKey, ExecutionStatus.Done, Option(0))
        } yield ()
        Await.result(createScatter, AkkaTimeout)
        newEntries.keys
      }
      case Success(v: WdlValue) => Failure(new Throwable("Scatter collection must evaluate to an array"))
      case Failure(ex) => Failure(ex)
    }
  }

  private def processRunnableCollector(collector: CollectorKey): Try[Iterable[ExecutionStoreKey]] = {
    persistStatus(collector, ExecutionStatus.Starting, None)
    val shards: Iterable[CallKey] = findShardEntries(collector) collect { case (k: CallKey, _) => k }

    generateCollectorOutput(collector, shards) match {
      case Failure(e) =>
        self ! CallFailed(collector, None, e.getMessage)
      case Success(outputs) =>
        log.info(s"Collection complete for Scattered Call ${collector.tag}.")
        self ! CallCompleted(collector, outputs, 0)
    }

    Success(Seq.empty[ExecutionStoreKey])
  }

  private def processRunnableCall(callKey: CallKey, status: ExecutionStatus): Try[Iterable[ExecutionStoreKey]] = {
    persistStatus(callKey, ExecutionStatus.Starting, None)
    val futureCallInputs = for {
      allInputs <- fetchLocallyQualifiedInputs(callKey)
    } yield allInputs
    Await.ready(futureCallInputs, AkkaTimeout)
    futureCallInputs.value.get match {
      case Success(callInputs) =>
        startActor(callKey, callInputs)
        Success(Seq.empty[ExecutionStoreKey])
      case Failure(ex) => Failure(ex)
    }
  }

  private def symbolsAsTable: Future[Traversable[Seq[String]]] = {
    def buildColumnData(entry: SymbolStoreEntry): Seq[String] = {
      // If the value of this symbol store entry is defined, that value truncated to MarkdownMaxColumnChars,
      // otherwise the empty string.
      def columnValue(entry: SymbolStoreEntry): String = {
        entry.wdlValue map { value =>
          s"(${value.wdlType.toWdlString}) ${value.valueString}".take(MarkdownMaxColumnChars)
        } getOrElse ""
      }
      // The Markdown columns as a Seq[String].
      Seq(
        entry.key.scope,
        entry.key.name,
        entry.key.index.map(_.toString).getOrElse(""),
        if (entry.key.input) "INPUT" else "OUTPUT",
        entry.wdlType.toWdlString,
        columnValue(entry)
      )
    }

    for {
      symbols <- globalDataAccess.getAllSymbolStoreEntries(workflow.id)
      columnData = symbols map buildColumnData
    } yield columnData
  }

  private def symbolsMarkdownTable: Future[Option[String]] = {
    val header = Seq("SCOPE", "NAME", "INDEX", "I/O", "TYPE", "VALUE")
    symbolsAsTable map {
      case rows if rows.isEmpty => None
      case rows => Option(TerminalUtil.mdTable(rows.toSeq, header))
    }
  }

  private def executionsAsTable: Future[Iterable[Seq[String]]] = {
    def buildColumnData(statusEntry: (ExecutionDatabaseKey, CallStatus)): Seq[String] = {
      val (k, v) = statusEntry
      Seq(k.fqn.toString, k.index.getOrElse("").toString, v.executionStatus.toString, v.returnCode.getOrElse("").toString)
    }

    for {
      executionStatuses <- globalDataAccess.getExecutionStatuses(workflow.id)
      columnData = executionStatuses map buildColumnData
    } yield columnData
  }

  private def executionsMarkdownTable: Future[Option[String]] = {
    val header = Seq("SCOPE", "INDEX", "STATUS", "RETURN CODE")
    executionsAsTable map {
      case rows if rows.isEmpty => None
      case rows => Option(TerminalUtil.mdTable(rows.toSeq, header))
    }
  }
}
