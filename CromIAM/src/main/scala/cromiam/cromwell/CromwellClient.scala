package cromiam.cromwell

import java.net.URL

import akka.actor.{ActorRef, ActorSystem}
import akka.event.LoggingAdapter
import akka.http.scaladsl.Http
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import akka.http.scaladsl.model.{HttpRequest, HttpResponse}
import akka.stream.ActorMaterializer
import com.softwaremill.sttp._
import cromiam.auth.{Collection, User}
import cromiam.cromwell.CromwellClient._
import cromiam.instrumentation.CromIamInstrumentation
import cromwell.api.{CromwellClient => CromwellApiClient}
import cromiam.server.status.StatusCheckedSubsystem
import cromwell.api.model.{WorkflowId, WorkflowLabels, WorkflowMetadata}
import spray.json._
import scala.concurrent.duration._
import scala.concurrent.{ExecutionContextExecutor, Future}
import scala.util.{Failure, Success}

/**
  * Provides a CromIAM specific handle for Cromwell communication
  *
  * FIXME: Look for ways to synch this up with the mothership
  */
class CromwellClient(scheme: String, interface: String, port: Int, log: LoggingAdapter, serviceRegistryActorRef: ActorRef)(implicit system: ActorSystem,
                                                                                        ece: ExecutionContextExecutor,
                                                                                        materializer: ActorMaterializer)
  extends SprayJsonSupport with DefaultJsonProtocol with StatusCheckedSubsystem with CromIamInstrumentation{

  val cromwellUrl = new URL(s"$scheme://$interface:$port")
  val cromwellApiVersion = "v1"

  override val statusUri = uri"$cromwellUrl/engine/$cromwellApiVersion/status"

  override val serviceRegistryActor: ActorRef = serviceRegistryActorRef

  val cromwellApiClient: CromwellApiClient = new CromwellApiClient(cromwellUrl, cromwellApiVersion)

  def collectionForWorkflow(workflowId: String, user: User): Future[Collection] = {
    import CromwellClient.EnhancedWorkflowLabels

    log.info("Requesting collection for " + workflowId + " for user " + user.userId + " from metadata")

    // Look up in Cromwell what the collection is for this workflow. If it doesn't exist, fail the Future
    cromwellApiClient.labels(WorkflowId.fromString(workflowId), headers = List(user.authorization)) flatMap {
      _.caasCollection match {
        case Some(c) => Future.successful(c)
        case None => Future.failed(new IllegalArgumentException(s"Workflow $workflowId has no associated collection"))
      }
    }
  }

  def forwardToCromwell(httpRequest: HttpRequest): Future[HttpResponse] = {
    val headers = httpRequest.headers.filterNot(header => header.name == TimeoutAccessHeader)
    val cromwellRequest = httpRequest
      .copy(uri = httpRequest.uri.withAuthority(interface, port).withScheme(scheme))
      .withHeaders(headers)
    Http().singleRequest(cromwellRequest)
  } recoverWith {
    case e => Future.failed(CromwellConnectionFailure(e))
  }

  /**
    * Retrieve the root workflow ID for a workflow ID. This is used in case the user is inquiring about a subworkflow.
    */
  def getRootWorkflow(workflowId: String, user: User, httpRequest: HttpRequest): Future[String] = {
    def metadataToRootWorkflowId(metadata: WorkflowMetadata): String = {
      import spray.json._
      /*
          Parse the JSON and then look for a top level field in the response object called rootWorkflowId, and return.
          If for some reason there's not a root workflow ID (not a subworkflow and/or this is an old workflow), just
          use the current workflow id.

          This is all called from inside the context of a Future, so exceptions will be properly caught.
      */
      metadata.value.parseJson.asJsObject.fields.get("rootWorkflowId").map(_.convertTo[String]).getOrElse(workflowId)
    }

    log.info("Looking up root workflow ID for " + workflowId + "for user " + user.userId + " from metadata")

    val startTimestamp = System.currentTimeMillis

    /*
      Grab the metadata from Cromwell filtered down to the rootWorkflowId. Then transform the response to get just the
      root workflow ID itself
     */
    val cromwellResponse = cromwellApiClient.metadata(WorkflowId.fromString(workflowId),
                              args = Option(Map("includeKey" -> List("rootWorkflowId"))),
                              headers = List(user.authorization)).map(metadataToRootWorkflowId)

    cromwellResponse.onComplete {
      case Success(_) => sendTimingApi(makeRequestPath(httpRequest, "success"), (System.currentTimeMillis - startTimestamp).millis, "root-workflow-id")
      case Failure(_) => sendTimingApi(makeRequestPath(httpRequest, "failed"), (System.currentTimeMillis - startTimestamp).millis, "root-workflow-id")
    }

    cromwellResponse
  }
}

object CromwellClient {
  // Header that Akka HTTP adds to every request on receive.
  // We get an warning in logs if we don't strip it out before sending the request to cromwell
  // HTTP header ‘Timeout-Access: <function1>’ is not allowed in requests
  // See: https://github.com/akka/akka-http/issues/64
  val TimeoutAccessHeader = "Timeout-Access"

  final case class CromwellConnectionFailure(f: Throwable) extends Exception(s"Unable to connect to Cromwell (${f.getMessage})", f)

  implicit class EnhancedWorkflowLabels(val wl: WorkflowLabels) extends AnyVal {

    import Collection.CollectionLabelName
    import Collection.collectionJsonReader

    def caasCollection: Option[Collection] = {
      wl.labels.fields.get(CollectionLabelName).map(_.convertTo[Collection])
    }
  }
}
