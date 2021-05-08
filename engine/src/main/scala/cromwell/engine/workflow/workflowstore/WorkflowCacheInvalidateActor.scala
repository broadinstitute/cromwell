package cromwell.engine.workflow.workflowstore

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.core._
import cromwell.engine.instrumentation.WorkflowInstrumentation
import cromwell.engine.workflow.WorkflowMetadataHelper
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCache
import cromwell.engine.workflow.workflowstore.WorkflowCacheInvalidateActor._
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor._
import cromwell.services.EngineServicesStore

import scala.concurrent.ExecutionContext
import scala.util.{Failure, Success}

final case class WorkflowCacheInvalidateActor(serviceRegistryActor: ActorRef) extends Actor
  with ActorLogging with WorkflowMetadataHelper with MonitoringCompanionHelper with WorkflowInstrumentation {
  implicit val ec: ExecutionContext = context.dispatcher
  private val callCache = new CallCache(EngineServicesStore.engineDatabaseInterface)

  val workflowCacheInvalidateReceive: Receive = {
    case cmd: InvalidateWorkflowCallCache =>
      addWork()
      val sndr = sender()

      callCache.invalidate(cmd.id) onComplete {
        case Success(()) =>
          log.info("workflow {} caches invalidated", cmd.id)
          sndr ! WorkflowCacheInvalidateSuccess
          removeWork()
        case Failure(throwable) =>
          log.error("Workflow {} cache invalidate failed", throwable)
          sndr ! WorkflowCacheInvalidateFailed(throwable)
          removeWork()
      }
  }

  override def receive: Receive = workflowCacheInvalidateReceive.orElse(monitoringReceive)
}

object WorkflowCacheInvalidateActor {
  def props(serviceRegistryActor: ActorRef): Props = {
    Props(WorkflowCacheInvalidateActor(serviceRegistryActor)).withDispatcher(ApiDispatcher)
  }

  sealed trait WorkflowCacheInvalidateActorResponse

  final case object WorkflowCacheInvalidateSuccess extends WorkflowCacheInvalidateActorResponse

  final case class WorkflowCacheInvalidateFailed(throwable: Throwable)
    extends WorkflowCacheInvalidateActorResponse
}
