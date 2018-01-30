package cromiam.cromwell

import java.net.URL

import akka.actor.ActorSystem
import akka.event.LoggingAdapter
import akka.http.scaladsl.Http
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import akka.http.scaladsl.model.{HttpRequest, HttpResponse}
import akka.stream.ActorMaterializer
import com.softwaremill.sttp._
import cromiam.auth.{Collection, User}
import cromiam.cromwell.CromwellClient._
import cromwell.api.{CromwellClient => CromwellApiClient}
import cromiam.server.status.StatusCheckedSubsystem
import cromwell.api.model.{WorkflowId, WorkflowLabels}
import spray.json._

import scala.concurrent.{ExecutionContextExecutor, Future}

/**
  * Provides a CromIAM specific handle for Cromwell communication
  *
  * FIXME: Look for ways to synch this up with the mothership
  */
class CromwellClient(scheme: String, interface: String, port: Int, log: LoggingAdapter)(implicit system: ActorSystem,
                                                                                        ece: ExecutionContextExecutor,
                                                                                        materializer: ActorMaterializer)
  extends SprayJsonSupport with DefaultJsonProtocol with StatusCheckedSubsystem {

  val cromwellUrl = new URL(s"$scheme://$interface:$port")
  val cromwellApiVersion = "v1"

  override val statusUri = uri"$cromwellUrl/engine/$cromwellApiVersion/status"

  /**
    FIXME:

    Creating a new client here every time sucks but currently headers on the CromwellClient are per-client and
    we need them to be per-request. A minor change to Cromwell, but beyond the skateboard we currently have
  */
  def collectionForWorkflow(workflowId: String, user: User): Future[Collection] = {
    import CromwellClient.EnhancedWorkflowLabels

    log.info("Requesting collection for " + workflowId + " for user " + user.userId + " from metadata")

    val client = new CromwellApiClient(cromwellUrl, cromwellApiVersion, Option(user.authorization.credentials))

    // Look up in Cromwell what the collection is for this workflow. If it doesn't exist, fail the Future
    client.labels(WorkflowId.fromString(workflowId)) flatMap {
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
      wl.labels.fields.get(CollectionLabelName) map { _.convertTo[Collection] }
    }
  }
}
