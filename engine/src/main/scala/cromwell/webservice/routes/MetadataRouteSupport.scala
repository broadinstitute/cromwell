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
import cromwell.core.{WorkflowId, path => _}
import cromwell.engine.instrumentation.HttpInstrumentation
import cromwell.server.CromwellShutdown
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.impl.builder.MetadataBuilderActor.{BuiltMetadataResponse, FailedMetadataResponse, MetadataBuilderActorResponse}
import cromwell.webservice.LabelsManagerActor
import cromwell.webservice.LabelsManagerActor._
import cromwell.webservice.routes.CromwellApiService.{InvalidWorkflowException, UnrecognizedWorkflowException, serviceShuttingDownResponse, validateWorkflowIdInMetadata, validateWorkflowIdInMetadataSummaries}
import cromwell.webservice.routes.MetadataRouteSupport._
import cromwell.webservice.WebServiceUtils.EnhancedThrowable
import cromwell.webservice.WorkflowJsonSupport._

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

            metadataLookup(possibleWorkflowId,
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
                case Valid(labels) =>
                  val response = validateWorkflowIdInMetadataSummaries(possibleWorkflowId, serviceRegistryActor) flatMap { id =>
                    val lma = actorRefFactory.actorOf(LabelsManagerActor.props(serviceRegistryActor).withDispatcher(ApiDispatcher))
                    lma.ask(LabelsAddition(LabelsData(id, labels))).mapTo[LabelsManagerActorResponse]
                  }
                  onComplete(response) {
                    case Success(r: BuiltLabelsManagerResponse) => complete(r.response)
                    case Success(e: FailedLabelsManagerResponse) => e.reason.failRequest(StatusCodes.InternalServerError)
                    case Failure(e: UnrecognizedWorkflowException) => e.failRequest(StatusCodes.NotFound)
                    case Failure(e: TimeoutException) => e.failRequest(StatusCodes.ServiceUnavailable)
                    case Failure(e) => e.errorRequest(StatusCodes.InternalServerError)

                  }
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
  def metadataLookup(possibleWorkflowId: String,
                     request: WorkflowId => MetadataReadAction,
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
                                  request: WorkflowId => MetadataReadAction,
                                  serviceRegistryActor: ActorRef)
                                 (implicit timeout: Timeout,
                                  ec: ExecutionContext): Future[MetadataBuilderActorResponse] = {
    validateWorkflowIdInMetadata(possibleWorkflowId, serviceRegistryActor) flatMap { w => serviceRegistryActor.ask(request(w)).mapTo[MetadataBuilderActorResponse] }
  }

  def completeMetadataBuilderResponse(response: Future[MetadataBuilderActorResponse]): Route = {
    onComplete(response) {
      case Success(r: BuiltMetadataResponse) => complete(r.responseJson)
      case Success(r: FailedMetadataResponse) => r.reason.errorRequest(StatusCodes.InternalServerError)
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
}
