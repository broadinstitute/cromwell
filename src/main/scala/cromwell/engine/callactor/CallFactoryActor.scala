//package cromwell.engine.callactor
//
//import akka.actor.{PoisonPill, ActorRef}
//import akka.event.Logging
//import akka.pattern.pipe
//import cromwell.engine._
//import cromwell.engine.backend.Backend
//import cromwell.engine.callactor.CallFactoryActor.{CallFactoryActorData, CallFactoryActorState}
//import cromwell.engine.db.DataAccess._
//import cromwell.engine.db.ExecutionDatabaseKey
//import cromwell.engine.db.slick.Execution
//import cromwell.engine.workflow.{ExecutionStoreKey, BackendCallKey}
//import cromwell.engine.workflow.WorkflowActor.{CallFailedToInitialize, InitialStartCall, UseCachedCall}
//import cromwell.logging.WorkflowLogger
//import wdl4s.Call
//import wdl4s.values.WdlValue
//
//import scala.concurrent.Future
//import scala.util.{Failure, Success, Try}
//
//object CallFactoryActor {
//  sealed trait CallFactoryActorState
//  case object Idle extends CallFactoryActorState
//  case object PersistingRuntimeAttributes extends CallFactoryActorState
//  case object WaitingForHash extends CallFactoryActorState
//  case object WaitingForReusableExecutions extends CallFactoryActorState
//  case object Done extends CallFactoryActorState
//
//  class CallFactoryActorData
//  sealed trait CallFactoryActorMessage
//  case object Start extends CallFactoryActorMessage
//  case object RuntimeAttributesPersisted extends CallFactoryActorMessage
//  case class AsyncFailure(failureContext: String, failure: Throwable) extends CallFactoryActorMessage
//  case class CallHash(hash: ExecutionHash)
//  case class ReusableExecutions(executions: Traversable[Execution])
//
//  implicit class EnhancedExecutionStoreKey(val key: ExecutionStoreKey) extends AnyVal {
//    def toDatabaseKey: ExecutionDatabaseKey = ExecutionDatabaseKey(key.scope.fullyQualifiedName, key.index, key.attempt)
//  }
//}
//
//
//class CallFactoryActor(workflow: WorkflowDescriptor, backend: Backend, callKey: BackendCallKey, callInputs: Map[String, WdlValue], replyTo: ActorRef) extends CromwellFSM[CallFactoryActorState, CallFactoryActorData] {
//  import CallFactoryActor._
//  val akkaLogger = Logging(context.system, classOf[CallFactoryActor])
//  implicit val logger: WorkflowLogger = WorkflowLogger("WorkflowActor", workflow, Option(akkaLogger))
//  val backendCall = backend.bindCall(workflow, callKey, callInputs, Option(AbortRegistrationFunction(registerAbortFunction)))
//  def registerAbortFunction(abortFunction: AbortFunction): Unit = {}
//
//  when(Idle) {
//    case Event(Start, _)=>
//      val validateAndPersistRA = for {
//        attrs <- Future.fromTry(Try(backendCall.runtimeAttributes))
//        globalDataAccess.setRuntimeAttributes(workflow.id, backendCall.key.toDatabaseKey, attrs.attributes)
//      } yield ()
//
//      chain(validateAndPersistRA map { _ => RuntimeAttributesPersisted }, PersistingRuntimeAttributes, s"Could not validate or persist runtime attributes for ${backendCall.key.tag}")
//  }
//
//  when(PersistingRuntimeAttributes) {
//    case Event(RuntimeAttributesPersisted, _) =>
//      if (backendCall.workflowDescriptor.readFromCache) {
//        val result = for {
//          hash <- backendCall.hash
//          executions <- globalDataAccess.getExecutionsWithResuableResultsByHash(hash.overallHash)
//        } yield ReusableExecutions(executions)
//
//        chain(result, WaitingForReusableExecutions, s"Failed to retrieve executions for reusable ${backendCall.key.tag}")
//      } else {
//        replyTo ! InitialStartCall(callKey, CallActor.Start)
//        self ! PoisonPill
//        goto(Done)
//      }
//  }
//
//  when(WaitingForReusableExecutions) {
//    case Event(ReusableExecutions(executions), _) =>
//      val result = globalDataAccess.getExecutionsWithResuableResultsByHash(hash.overallHash) map ReusableExecutions.apply recover {
//        case e => ActorFailure(s"Call Caching: Failed to look up executions that matched hash '$hash'. Falling back to normal execution", e) }
//      }
//        val reusableExecutions = globalDataAccess.getExecutionsWithResuableResultsByHash(hash.overallHash) onComplete {
//          case Success(executions) if executions.nonEmpty => startCachedCall(executions.head)
//          case Success(_) =>
//            log.info(s"Call Caching: cache miss")
//            self ! InitialStartCall(callKey, CallActor.Start)
//          case Failure(ex) =>
//            log.error(s"Call Caching: Failed to look up executions that matched hash '$hash'. Falling back to normal execution", ex)
//            self ! InitialStartCall(callKey, CallActor.Start)
//        }
//  }
//
//    def loadCachedBackendCallAndMessage(descriptor: WorkflowDescriptor, cachedExecution: Execution) = {
//      descriptor.namespace.resolve(cachedExecution.callFqn) match {
//        case Some(c: Call) =>
//          val cachedCall = backend.bindCall(
//            descriptor,
//            BackendCallKey(c, cachedExecution.index.toIndex, cachedExecution.attempt),
//            callInputs,
//            Option(AbortRegistrationFunction(registerAbortFunction))
//          )
//          log.info(s"Call Caching: Cache hit. Using UUID(${cachedCall.workflowDescriptor.shortId}):${cachedCall.key.tag} as results for UUID(${backendCall.workflowDescriptor.shortId}):${backendCall.key.tag}")
//          self ! UseCachedCall(callKey, CallActor.UseCachedCall(cachedCall, backendCall))
//        case _ =>
//          log.error(s"Call Caching: error when resolving '${cachedExecution.callFqn}' in workflow with execution ID ${cachedExecution.workflowExecutionId}: falling back to normal execution")
//          self ! InitialStartCall(callKey, CallActor.Start)
//      }
//    }
//
//    /* Tries to use the cached Execution to send a UseCachedCall message.  If anything fails, send an InitialStartCall message */
//    def loadCachedCallOrInitiateCall(cachedDescriptor: Try[WorkflowDescriptor], cachedExecution: Execution) = cachedDescriptor match {
//      case Success(descriptor) => loadCachedBackendCallAndMessage(descriptor, cachedExecution)
//      case Failure(ex) =>
//        log.error(s"Call Caching: error when loading workflow with execution ID ${cachedExecution.workflowExecutionId}: falling back to normal execution", ex)
//        self ! InitialStartCall(callKey, CallActor.Start)
//    }
//
//    def startCachedCall(cachedExecution: Execution) = {
//      globalDataAccess.getWorkflow(cachedExecution.workflowExecutionId) onComplete { cachedDescriptor =>
//        loadCachedCallOrInitiateCall(cachedDescriptor, cachedExecution)
//      }
//    }
//
//    def checkCacheAndStartCall = {
//      log.debug(s"Call caching 'readFromCache' is turned on. Checking cache before starting call")
//      val result = backendCall.hash map CallHash.apply recover { case e => ActorFailure(s"Failed to calculate hash for call '${backendCall.key.tag}'.", e) }
//      pipe(result) to self
//
//      //        globalDataAccess.getExecutionsWithResuableResultsByHash(hash.overallHash) onComplete {
//      //          case Success(executions) if executions.nonEmpty => startCachedCall(executions.head)
//      //          case Success(_) =>
//      //            log.info(s"Call Caching: cache miss")
//      //            self ! InitialStartCall(callKey, CallActor.Start)
//      //          case Failure(ex) =>
//      //            log.error(s"Call Caching: Failed to look up executions that matched hash '$hash'. Falling back to normal execution", ex)
//      //            self ! InitialStartCall(callKey, CallActor.Start)
//      //        }
//    }
//
//}
