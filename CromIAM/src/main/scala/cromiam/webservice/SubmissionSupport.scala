package cromiam.webservice

import akka.http.scaladsl.model.{HttpEntity, _}
import akka.http.scaladsl.server._
import cromiam.cromwell.CromwellClient
import cromiam.sam.SamClient
import cromiam.sam.SamClient._
import akka.http.scaladsl.server.Directives._
import akka.stream.ActorMaterializer
import akka.util.ByteString
import SubmissionSupport._
import cromiam.auth.{Collection, User}
import akka.event.LoggingAdapter
import spray.json.{JsObject, JsString, JsValue}
import Collection.{CollectionLabelName, LabelsKey, validateLabels}

import scala.concurrent.{ExecutionContextExecutor, Future}

trait SubmissionSupport extends RequestSupport {
  val cromwellClient: CromwellClient
  val samClient: SamClient

  val log: LoggingAdapter

  implicit def executor: ExecutionContextExecutor
  implicit val materializer: ActorMaterializer

  // FIXME - getting pathPrefix to shrink this keeps hosing up, there's gotta be some way to do this
  def submitRoute: Route = (path("api" / "workflows" / Segment) | path("api" / "workflows" / Segment / "batch")) { _ =>
    post {
      extractUserAndRequest { (user, request) =>
        log.info("Received submission request from user " + user.userId)
        authorizeAsync(samClient.isSubmitWhitelisted(user)) {
          extractSubmission(user) { submission =>
            complete {
              forwardSubmissionToCromwell(user, submission.collection, request.withEntity(submission.entity))
            }
          }
        }
      }
    }
  }

  private def forwardSubmissionToCromwell(user: User,
                                          collection: Collection,
                                          submissionRequest: HttpRequest): Future[HttpResponse] = {
    log.info("Forwarding submission request for " + user.userId + " with collection " + collection.name + " to Cromwell")

    def registerWithSam(collection: Collection): Future[Unit] = {
      samClient.requestSubmission(user, collection) recoverWith {
        case SamDenialException => Future.failed(SamDenialException)
        case e => Future.failed(SamConnectionFailure("new workflow registration", e))
      }
    }

    (for {
      resp <- cromwellClient.forwardToCromwell(submissionRequest)
      _ <- registerWithSam(collection)
    } yield resp) recover {
      case SamDenialException => SamDenialResponse
    }
  }
}

object SubmissionSupport {
  def extractCollection(user: User): Directive1[Collection] = {
    formField(CollectionNameKey.?) map { maybeCollectionName =>
      maybeCollectionName.map(Collection(_)).getOrElse(Collection.forUser(user))
    }
  }

  def extractSubmission(user: User): Directive1[WorkflowSubmission] = {
    (
      extractCollection(user) &
      formFields((
        WorkflowSourceKey.?,
        WorkflowUrlKey.?,
        WorkflowTypeKey.?,
        WorkflowTypeVersionKey.?,
        WorkflowInputsKey.?,
        WorkflowOptionsKey.?,
        WorkflowOnHoldKey.as[Boolean].?,
        WorkflowDependenciesKey.as[ByteString].?)) &
      extractLabels &
      extractInputAux
    ).as(WorkflowSubmission)
  }

  def extractLabels: Directive1[Option[Map[String, JsValue]]] = {
    formField(LabelsKey.?) flatMap validateLabels
  }

  def extractInputAux: Directive1[Map[String, String]] = {
    formFieldMap.map(_.filterKeys(_.startsWith(WorkflowInputsAuxPrefix)))
  }

  // FIXME: Much like CromwellClient see if there are ways of unifying this a bit w/ the mothership
  final case class WorkflowSubmission(collection: Collection,
                                      workflowSource: Option[String],
                                      workflowUrl: Option[String],
                                      workflowType: Option[String],
                                      workflowTypeVersion: Option[String],
                                      workflowInputs: Option[String],
                                      workflowOptions: Option[String],
                                      workflowOnHold: Option[Boolean],
                                      workflowDependencies: Option[ByteString],
                                      origLabels: Option[Map[String, JsValue]],
                                      workflowInputsAux: Map[String, String]) {
    // For auto-validation, if origLabels defined, can't have CaaS collection label set. Was checked previously, but ...
    require(origLabels.forall(!_.keySet.contains(CollectionLabelName)))

    // Inject the collection name into the labels and convert to a String
    private val collectionLabels = Map(CollectionLabelName -> JsString(collection.name))
    private val labels: String = JsObject(origLabels.map(o => o ++ collectionLabels).getOrElse(collectionLabels)).toString

    val entity: MessageEntity = {
      val sourcePart = workflowSource map { s => Multipart.FormData.BodyPart(WorkflowSourceKey, HttpEntity(MediaTypes.`application/json`, s)) }
      val urlPart = workflowUrl map { u => Multipart.FormData.BodyPart(WorkflowUrlKey, HttpEntity(MediaTypes.`application/json`, u))}
      val typePart = workflowType map { t => Multipart.FormData.BodyPart(WorkflowTypeKey, HttpEntity(MediaTypes.`application/json`, t)) }
      val typeVersionPart = workflowTypeVersion map { v => Multipart.FormData.BodyPart(WorkflowTypeVersionKey, HttpEntity(MediaTypes.`application/json`, v)) }
      val inputsPart = workflowInputs map { i => Multipart.FormData.BodyPart(WorkflowInputsKey, HttpEntity(MediaTypes.`application/json`, i)) }
      val optionsPart = workflowOptions map { o => Multipart.FormData.BodyPart(WorkflowOptionsKey, HttpEntity(MediaTypes.`application/json`, o)) }
      val importsPart = workflowDependencies map { d => Multipart.FormData.BodyPart(WorkflowDependenciesKey, HttpEntity(MediaTypes.`application/octet-stream`, d)) }
      val onHoldPart = workflowOnHold map { h => Multipart.FormData.BodyPart(WorkflowOnHoldKey, HttpEntity(h.toString)) }
      val labelsPart = Multipart.FormData.BodyPart(LabelsKey, HttpEntity(MediaTypes.`application/json`, labels))
      val parts = List(sourcePart, urlPart, typePart, typeVersionPart, inputsPart, optionsPart, importsPart, onHoldPart, Option(labelsPart)).flatten ++ auxParts

      Multipart.FormData(parts: _*).toEntity()
    }

    private def auxParts = {
      workflowInputsAux map { case (k, v) => Multipart.FormData.BodyPart(k, HttpEntity(MediaTypes.`application/json`, v)) }
    }
  }

  // FIXME: Unify these w/ Cromwell.PartialWorkflowSources (via common?)
  val CollectionNameKey = "collectionName"
  val WorkflowSourceKey = "workflowSource"
  val WorkflowUrlKey = "workflowUrl"
  val WorkflowTypeKey = "workflowType"
  val WorkflowTypeVersionKey = "workflowTypeVersion"
  val WorkflowInputsKey = "workflowInputs"
  val WorkflowInputsAuxPrefix = "workflowInputs_"
  val WorkflowOptionsKey = "workflowOptions"
  val WorkflowDependenciesKey = "workflowDependencies"
  val WorkflowOnHoldKey = "workflowOnHold"
}
