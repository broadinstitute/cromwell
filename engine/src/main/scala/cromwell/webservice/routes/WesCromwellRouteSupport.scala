package cromwell.webservice.routes


import akka.actor.{ActorRef, ActorRefFactory}
import akka.http.javadsl.server.Directives.handleExceptions
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model.{Multipart, StatusCodes}
import akka.http.scaladsl.server.Directives.onComplete
import akka.http.scaladsl.server.{ExceptionHandler, Route}
import akka.pattern.{AskTimeoutException, ask}
import akka.stream.ActorMaterializer
import akka.util.Timeout
import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import cromwell.core.abort.SuccessfulAbortResponse
import cromwell.core.{WorkflowId, WorkflowOnHold, WorkflowState, WorkflowSubmitted, path => _}
import cromwell.engine.workflow.workflowstore.WorkflowStoreSubmitActor.{WorkflowStoreSubmitActorResponse, WorkflowSubmitFailed}
import cromwell.engine.workflow.workflowstore.{WorkflowStoreActor, WorkflowStoreSubmitActor}
import cromwell.server.CromwellShutdown
import cromwell.webservice.WebServiceUtils.EnhancedThrowable
import cromwell.webservice.WorkflowJsonSupport._
import cromwell.webservice.routes.CromwellApiService.standardAbortSuccessHandler
import cromwell.webservice.{PartialWorkflowSources, WebServiceUtils, WorkflowSubmitResponse}
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration.FiniteDuration
import scala.concurrent.{ExecutionContext, TimeoutException}
import scala.util.{Failure, Success}

trait WesCromwellRouteSupport extends WebServiceUtils {

  val workflowStoreActor: ActorRef


  implicit val duration = ConfigFactory.load().as[FiniteDuration]("akka.http.server.request-timeout")
  implicit val timeout: Timeout = duration

  implicit def actorRefFactory: ActorRefFactory

  implicit val materializer: ActorMaterializer
  implicit val ec: ExecutionContext

  def toResponse(workflowId: WorkflowId, workflowState: WorkflowState): WorkflowSubmitResponse = {
    WorkflowSubmitResponse(workflowId.toString, workflowState.toString)
  }

  def submitRequest(formData: Multipart.FormData, isSingleSubmission: Boolean): Route = {

    def getWorkflowState(workflowOnHold: Boolean): WorkflowState = {
      if (workflowOnHold)
        WorkflowOnHold
      else WorkflowSubmitted
    }

    def sendToWorkflowStore(command: WorkflowStoreActor.WorkflowStoreActorSubmitCommand, warnings: Seq[String], workflowState: WorkflowState): Route = {
      // NOTE: Do not blindly copy the akka-http -to- ask-actor pattern below without knowing the pros and cons.
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
        case Failure(_: AskTimeoutException) if CromwellShutdown.shutdownInProgress() => CromwellApiService.serviceShuttingDownResponse
        case Failure(e: TimeoutException) => e.failRequest(StatusCodes.ServiceUnavailable)
        case Failure(e) => e.failRequest(StatusCodes.InternalServerError, warnings)
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
