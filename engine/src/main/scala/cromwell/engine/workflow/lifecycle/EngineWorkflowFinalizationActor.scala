package cromwell.engine.workflow.lifecycle

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.event.LoggingReceive
import cromwell.backend.BackendLifecycleActor.BackendWorkflowLifecycleActorResponse
import cromwell.backend.BackendWorkflowFinalizationActor.{FinalizationFailed, FinalizationResponse, Finalize}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

trait EngineWorkflowFinalizationActor extends Actor with ActorLogging {
  def receive: Receive = LoggingReceive {
    case Finalize => performActionThenRespond(afterAll()(context.dispatcher), FinalizationFailed)(context.dispatcher)
  }

  protected def performActionThenRespond(operation: => Future[BackendWorkflowLifecycleActorResponse],
                                         onFailure: (Throwable) => BackendWorkflowLifecycleActorResponse)
                                        (implicit ec: ExecutionContext) = {
    val respondTo: ActorRef = sender
    operation onComplete {
      case Success(r) => respondTo ! r
      case Failure(t) => respondTo ! onFailure(t)
    }
  }

  /**
    * Happens after everything else runs
    */
  def afterAll()(implicit ec: ExecutionContext): Future[FinalizationResponse]
}
