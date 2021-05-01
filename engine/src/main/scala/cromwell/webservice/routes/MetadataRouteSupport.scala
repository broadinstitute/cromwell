package cromwell.webservice.routes

import akka.actor.{ActorRef, ActorRefFactory}
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.marshalling.ToResponseMarshallable
import akka.http.scaladsl.model.StatusCodes
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.Route
import akka.pattern.{AskTimeoutException, ask}
import akka.util.Timeout
import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.core.labels.Labels
import cromwell.core.{WorkflowId, WorkflowMetadataKeys, path => _}
import cromwell.engine.instrumentation.HttpInstrumentation
import cromwell.server.CromwellShutdown
import cromwell.services._
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataService._
import cromwell.webservice.LabelsManagerActor
import cromwell.webservice.LabelsManagerActor._
import cromwell.webservice.WebServiceUtils.EnhancedThrowable
import cromwell.webservice.WorkflowJsonSupport._
import cromwell.webservice.routes.CromwellApiService._
import cromwell.webservice.routes.MetadataRouteSupport._
import cromwell.webservice.WebServiceUtils._
import spray.json.{JsObject, JsString}

import scala.concurrent.{ExecutionContext, Future, TimeoutException}
import scala.util.{Failure, Success}


trait MetadataRouteSupport extends HttpInstrumentation {
  implicit def actorRefFactory: ActorRefFactory
  implicit val ec: ExecutionContext

  val serviceRegistryActor: ActorRef

  implicit val timeout: Timeout

  val metadataRoutes = concat(
    path("workflows" / Segment / Segment / "status") { (_, possibleWorkflowId) =>
      get {
        instrumentRequest {
          metadataLookup(possibleWorkflowId, (w: WorkflowId) => GetStatus(w), serviceRegistryActor)
        }
      }
    },
    path("workflows" / Segment / Segment / "outputs") { (_, possibleWorkflowId) =>
      get {
        instrumentRequest {
          metadataLookup(possibleWorkflowId, (w: WorkflowId) => WorkflowOutputs(w), serviceRegistryActor)
        }
      }
    },
    path("workflows" / Segment / Segment / "logs") { (_, possibleWorkflowId) =>
      get {
        instrumentRequest {
          metadataLookup(possibleWorkflowId, (w: WorkflowId) => GetLogs(w), serviceRegistryActor)
        }
      }
    },
    encodeResponse {
      path("workflows" / Segment / Segment / "metadata") { (_, possibleWorkflowId) =>
        instrumentRequest {
          parameters(('includeKey.*, 'excludeKey.*, 'expandSubWorkflows.as[Boolean].?)) { (includeKeys, excludeKeys, expandSubWorkflowsOption) =>
            val includeKeysOption = NonEmptyList.fromList(includeKeys.toList)
            val excludeKeysOption = NonEmptyList.fromList(excludeKeys.toList)
            val expandSubWorkflows = expandSubWorkflowsOption.getOrElse(false)

            metadataLookup(
              possibleWorkflowId,
              (w: WorkflowId) => GetSingleWorkflowMetadataAction(w, includeKeysOption, excludeKeysOption, expandSubWorkflows),
              serviceRegistryActor)
          }
        }
      }
    },
    path("workflows" / Segment / Segment / "labels") { (_, possibleWorkflowId) =>
      concat(
        get {
          instrumentRequest {
            metadataLookup(possibleWorkflowId, (w: WorkflowId) => GetLabels(w), serviceRegistryActor)
          }
        },
        patch {
          entity(as[Map[String, String]]) { parameterMap =>
            instrumentRequest {
              Labels.validateMapOfLabels(parameterMap) match {
                case Valid(labels) => patchLabelsRequest(possibleWorkflowId, labels, serviceRegistryActor, actorRefFactory)
                case Invalid(e) =>
                  val iae = new IllegalArgumentException(e.toList.mkString(","))
                  iae.failRequest(StatusCodes.BadRequest)
              }
            }
          }
        }
      )
    },
    path("workflows" / Segment / "query") { _ =>
      get {
        instrumentRequest {
          parameterSeq { parameters =>
            queryMetadata(parameters, serviceRegistryActor)
          }
        }
      } ~
        post {
          instrumentRequest {
            entity(as[Seq[Map[String, String]]]) { parameterMap =>
              queryMetadata(parameterMap.flatMap(_.toSeq), serviceRegistryActor)
            }
          }
        }
    }
  )
}

object MetadataRouteSupport {

  private def processMetadataArchivedResponse(workflowId: WorkflowId,
                                              archiveStatus: MetadataArchiveStatus,
                                              additionalMsg: String = ""): JsObject = {
    val baseMessage = "Cromwell has archived this workflow's metadata according to the lifecycle policy."
    val additionalDetails = if (archiveStatus == MetadataArchiveStatus.ArchivedAndDeleted)
      " It is available in the archive bucket, or via a support request in the case of a managed instance." + additionalMsg
    else additionalMsg

    JsObject(Map(
      WorkflowMetadataKeys.Id -> JsString(workflowId.toString),
      WorkflowMetadataKeys.MetadataArchiveStatus -> JsString(archiveStatus.toString),
      WorkflowMetadataKeys.Message -> JsString(baseMessage + additionalDetails)
    ))
  }

  def metadataLookup(possibleWorkflowId: String,
                     request: WorkflowId => BuildMetadataJsonAction,
                     serviceRegistryActor: ActorRef)
                    (implicit timeout: Timeout,
                     ec: ExecutionContext): Route = {
    completeMetadataBuilderResponse(metadataBuilderActorRequest(possibleWorkflowId, request, serviceRegistryActor))
  }

  def queryMetadata(parameters: Seq[(String, String)],
                    serviceRegistryActor: ActorRef)(implicit timeout: Timeout): Route = {
    completeMetadataQueryResponse(metadataQueryRequest(parameters, serviceRegistryActor))
  }

  def metadataBuilderActorRequest(possibleWorkflowId: String,
                                  request: WorkflowId => BuildMetadataJsonAction,
                                  serviceRegistryActor: ActorRef)
                                 (implicit timeout: Timeout,
                                  ec: ExecutionContext): Future[MetadataJsonResponse] = {

    def checkIfMetadataDeletedAndRespond(id: WorkflowId,
                                         metadataRequest: BuildWorkflowMetadataJsonWithOverridableSourceAction): Future[MetadataJsonResponse] = {
      serviceRegistryActor.ask(FetchWorkflowMetadataArchiveStatus(id)).mapTo[FetchWorkflowArchiveStatusResponse] flatMap {
        case WorkflowMetadataArchivedStatus(archiveStatus) =>
          if (archiveStatus.isDeleted) Future.successful(SuccessfulMetadataJsonResponse(metadataRequest, processMetadataArchivedResponse(id, archiveStatus)))
          else serviceRegistryActor.ask(request(id)).mapTo[MetadataJsonResponse]
        case FailedToGetArchiveStatus(e) => Future.failed(e)
      }
    }

    validateWorkflowIdInMetadata(possibleWorkflowId, serviceRegistryActor) flatMap { id =>
      /*
        for requests made to one of /metadata, /logs or /outputs endpoints, perform an additional check to see
        if metadata for the workflow has been archived and deleted or not (as they interact with metadata table)
      */
      request(id) match {
        case m: BuildWorkflowMetadataJsonWithOverridableSourceAction => checkIfMetadataDeletedAndRespond(id, m)
        case _ => serviceRegistryActor.ask(request(id)).mapTo[MetadataJsonResponse]
      }
    }
  }

  def completeMetadataBuilderResponse(response: Future[MetadataJsonResponse]): Route = {
    onComplete(response) {
      case Success(r: SuccessfulMetadataJsonResponse) => complete(r.responseJson)
      case Success(r: FailedMetadataJsonResponse) => r.reason.errorRequest(StatusCodes.InternalServerError)
      case Failure(_: AskTimeoutException) if CromwellShutdown.shutdownInProgress() => serviceShuttingDownResponse
      case Failure(e: UnrecognizedWorkflowException) => e.failRequest(StatusCodes.NotFound)
      case Failure(e: InvalidWorkflowException) => e.failRequest(StatusCodes.BadRequest)
      case Failure(e: TimeoutException) => e.failRequest(StatusCodes.ServiceUnavailable)
      case Failure(e) => e.errorRequest(StatusCodes.InternalServerError)
    }
  }

  def metadataQueryRequest(parameters: Seq[(String, String)],
                           serviceRegistryActor: ActorRef)(implicit timeout: Timeout): Future[MetadataQueryResponse] = {
    serviceRegistryActor.ask(QueryForWorkflowsMatchingParameters(parameters)).mapTo[MetadataQueryResponse]
  }

  def completeMetadataQueryResponse(response: Future[MetadataQueryResponse]): Route = {
    import cromwell.webservice.WorkflowJsonSupport.workflowQueryResponse

    onComplete(response) {
      case Success(w: WorkflowQuerySuccess) => complete(ToResponseMarshallable(w.response))
      case Success(w: WorkflowQueryFailure) => w.reason.failRequest(StatusCodes.BadRequest)
      case Failure(_: AskTimeoutException) if CromwellShutdown.shutdownInProgress() => serviceShuttingDownResponse
      case Failure(e: TimeoutException) => e.failRequest(StatusCodes.ServiceUnavailable)
      case Failure(e) => e.errorRequest(StatusCodes.InternalServerError)
    }
  }

  def completePatchLabelsResponse(response: Future[LabelsManagerActorResponse]): Route = {
    onComplete(response) {
      case Success(r: BuiltLabelsManagerResponse) => complete(r.response)
      case Success(r: WorkflowArchivedLabelsManagerResponse) => completeResponse(StatusCodes.BadRequest, r.response, Seq.empty)
      case Success(e: FailedLabelsManagerResponse) => e.reason.failRequest(StatusCodes.InternalServerError)
      case Failure(e: UnrecognizedWorkflowException) => e.failRequest(StatusCodes.NotFound)
      case Failure(e: TimeoutException) => e.failRequest(StatusCodes.ServiceUnavailable)
      case Failure(e) => e.errorRequest(StatusCodes.InternalServerError)
    }
  }

  def patchLabelsRequest(possibleWorkflowId: String,
                         labels: Labels,
                         serviceRegistryActor: ActorRef,
                         actorRefFactory: ActorRefFactory)
                        (implicit timeout: Timeout, ec: ExecutionContext): Route = {

    def checkIfMetadataArchivedAndRespond(id: WorkflowId, archiveStatusResponse: FetchWorkflowArchiveStatusResponse): Future[LabelsManagerActorResponse] = {
      archiveStatusResponse match {
        case WorkflowMetadataArchivedStatus(archiveStatus) =>
          if (archiveStatus.isArchived) {
            val message = " As a result, new labels can't be added or existing labels can't be updated for this workflow."
            Future.successful(WorkflowArchivedLabelsManagerResponse(processMetadataArchivedResponse(id, archiveStatus, message)))
          }
          else {
            val lma = actorRefFactory.actorOf(LabelsManagerActor.props(serviceRegistryActor).withDispatcher(ApiDispatcher))
            lma.ask(LabelsAddition(LabelsData(id, labels))).mapTo[LabelsManagerActorResponse]
          }
        case FailedToGetArchiveStatus(e) => Future.failed(e)
      }
    }

    val response = for {
      id <- validateWorkflowIdInMetadataSummaries(possibleWorkflowId, serviceRegistryActor)
      archiveStatusResponse <- serviceRegistryActor.ask(FetchWorkflowMetadataArchiveStatus(id)).mapTo[FetchWorkflowArchiveStatusResponse]
      response <- checkIfMetadataArchivedAndRespond(id, archiveStatusResponse)
    } yield response

    completePatchLabelsResponse(response)
  }
}
