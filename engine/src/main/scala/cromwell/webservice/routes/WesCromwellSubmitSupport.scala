package cromwell.webservice.routes

import akka.actor.ActorRef
import akka.http.scaladsl.model.{Multipart, StatusCodes}
import akka.http.scaladsl.server.Directives.{handleExceptions, onComplete}
import akka.http.scaladsl.server.{ExceptionHandler, Route}
import akka.pattern.ask
import akka.util.Timeout
import cats.data.NonEmptyList
import net.ceedubs.ficus.Ficus._
import com.typesafe.config.ConfigFactory
import cromwell.core.{WorkflowOnHold, WorkflowState, WorkflowSubmitted}
import cromwell.engine.workflow.workflowstore.{WorkflowStoreActor, WorkflowStoreSubmitActor}
import cromwell.webservice.{PartialWorkflowSources, WebServiceUtils}
import cromwell.webservice.WebServiceUtils.EnhancedThrowable

import scala.concurrent.TimeoutException
import scala.concurrent.duration.FiniteDuration
import scala.util.{Failure, Success}

trait WesCromwellRouteSupport extends WebServiceUtils {

  val workflowStoreActor: ActorRef

  implicit val duration = ConfigFactory.load().as[FiniteDuration]("akka.http.server.request-timeout")
  implicit val timeout: Timeout = duration

  def submitRequest(formData: Multipart.FormData,
                    isSingleSubmission: Boolean,
                    //successHandler: PartialFunction[WorkflowStoreSubmitActorResponse, Route] = standardSuccessHandler,
                    //errorHandler: PartialFunction[Throwable, Route] = standardErrorHandler
                   ): Route = {

    def getWorkflowState(workflowOnHold: Boolean): WorkflowState = {
      if (workflowOnHold)
        WorkflowOnHold
      else WorkflowSubmitted
    }

    def sendToWorkflowStore(command: WorkflowStoreActor.WorkflowStoreActorSubmitCommand, warnings: Seq[String], workflowState: WorkflowState): Route = {
      //sendToWorkflowStore
      // NOTE: Do not blindly copy the akka-http -to- ask-actor pattern below without knowing the pros and cons.
      handleExceptions(ExceptionHandler(errorHandler)) {
        onComplete(workflowStoreActor.ask(command).mapTo[WorkflowStoreSubmitActor.WorkflowStoreSubmitActorResponse]) {
          case Success(w) =>
            w match {
              case WorkflowStoreSubmitActor.WorkflowSubmittedToStore(workflowId, _) =>
                completeResponse(StatusCodes.Created, toResponse(workflowId, workflowState), warnings)
              case WorkflowStoreSubmitActor.WorkflowsBatchSubmittedToStore(workflowIds, _) =>
                completeResponse(StatusCodes.Created, workflowIds.toList.map(toResponse(_, workflowState)), warnings)
              case WorkflowStoreSubmitActor.WorkflowSubmitFailed(throwable) =>
                throwable.failRequest(StatusCodes.BadRequest, warnings)
            }
          case Failure(e) => throw e
        }
      }
    }

    onComplete(materializeFormData(formData)) {
      case Success(data) =>
        PartialWorkflowSources.fromSubmitRoute(data, allowNoInputs = isSingleSubmission) match {
          case Success(workflowSourceFiles) if isSingleSubmission && workflowSourceFiles.size == 1 =>
            val warnings = workflowSourceFiles.flatMap(_.warnings)
            sendToWorkflowStore(WorkflowStoreActor.SubmitWorkflow(workflowSourceFiles.head), warnings, getWorkflowState(workflowSourceFiles.head.workflowOnHold))
          // Catches the case where someone has gone through the single submission endpoint w/ more than one workflow
          case Success(workflowSourceFiles) if isSingleSubmission =>
            val warnings = workflowSourceFiles.flatMap(_.warnings)
            val e = new IllegalArgumentException("To submit more than one workflow at a time, use the batch endpoint.")
            e.failRequest(StatusCodes.BadRequest, warnings)
          case Success(workflowSourceFiles) =>
            val warnings = workflowSourceFiles.flatMap(_.warnings)
            sendToWorkflowStore(
              WorkflowStoreActor.BatchSubmitWorkflows(NonEmptyList.fromListUnsafe(workflowSourceFiles.toList)),
              warnings, getWorkflowState(workflowSourceFiles.head.workflowOnHold))
          case Failure(t) => t.failRequest(StatusCodes.BadRequest)
        }
      case Failure(e: TimeoutException) => e.failRequest(StatusCodes.ServiceUnavailable)
      case Failure(e) => e.failRequest(StatusCodes.InternalServerError)
    }
  }
}

}
