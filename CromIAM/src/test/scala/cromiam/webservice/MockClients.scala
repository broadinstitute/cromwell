package cromiam.webservice

import akka.actor.{ActorRef, ActorSystem}
import akka.event.NoLogging
import akka.http.scaladsl.model.{HttpRequest, HttpResponse, StatusCodes}
import akka.http.scaladsl.model.StatusCodes._
import akka.stream.ActorMaterializer
import cromiam.auth.{Collection, User}
import cromiam.cromwell.CromwellClient
import cromiam.sam.SamClient
import cromiam.sam.SamClient.{CollectionAuthorizationRequest, SamDenialException, SamRegisterCollectionException}

import scala.concurrent.{ExecutionContextExecutor, Future}

class MockCromwellClient()(implicit system: ActorSystem,
                                 ece: ExecutionContextExecutor,
                                 materializer: ActorMaterializer)
  extends CromwellClient("http", "bar", 1, NoLogging, ActorRef.noSender)(system, ece, materializer) {
  val version = "v1"

  val unauthorizedUserCollectionStr: String = "987654321"
  val authorizedUserCollectionStr: String = "123456789"

  val subworkflowId = "58114f5c-f439-4488-8d73-092273cf92d9"
  val rootWorkflowIdWithCollection = "998d137e-7213-44ac-8f6f-24e6e23adaa5"
  val workflowIdWithoutCollection = "472234d-7213-44ac-8f6f-24e6e23adaa5"
  val anotherRootWorkflowIdWithCollection = "12b9a7c1-84b0-451e-96d7-a389315a2fa9"

  val userCollection = Collection("foo")

  def checkHeaderAuthorization(httpRequest: HttpRequest): Future[HttpResponse] = {
    val userIdHeader = httpRequest.headers.find(header => header.name.equalsIgnoreCase("OIDC_CLAIM_user_id"))

    userIdHeader match {
      case Some(header) =>
        header.value match {
          case s if s.equalsIgnoreCase(authorizedUserCollectionStr) => Future.successful(HttpResponse(status = OK, entity = "Response from Cromwell"))
          case s => Future.failed(new Exception(s"This is unexpected! Cromwell should not receive request from unauthorized user! OIDC_CLAIM_user_id: $s is not authorized."))
        }
      case None => Future.failed(new Exception("This is unexpected! No OIDC_CLAIM_user_id provided for authorization!"))
    }
  }

  override def forwardToCromwell(httpRequest: HttpRequest): Future[HttpResponse] = {
    val versionRoutePath = s"/engine/$version/version"

    httpRequest.uri.path.toString match {
      //version endpoint doesn't require authentication
      case `versionRoutePath` => Future.successful(HttpResponse(status = OK, entity = "Response from Cromwell"))
      case _ => checkHeaderAuthorization(httpRequest)
    }
  }

  override def getRootWorkflow(workflowId: String, user: User, cromIamRequest: HttpRequest): Future[String] = {
    workflowId match {
      case `subworkflowId` | `rootWorkflowIdWithCollection` => Future.successful(rootWorkflowIdWithCollection)
      case  `anotherRootWorkflowIdWithCollection` => Future.successful(anotherRootWorkflowIdWithCollection)
      case _ => Future.successful(workflowIdWithoutCollection)
    }
  }

  override def collectionForWorkflow(workflowId: String, user: User, cromIamRequest: HttpRequest): Future[Collection] = {
    workflowId match {
      case `rootWorkflowIdWithCollection` | `anotherRootWorkflowIdWithCollection` => Future.successful(userCollection)
      case _ => Future.failed(new IllegalArgumentException(s"Workflow $workflowId has no associated collection"))
    }
  }
}

class MockSamClient(checkSubmitWhitelist: Boolean = true)
                   (implicit system: ActorSystem,
                      ece: ExecutionContextExecutor,
                      materializer: ActorMaterializer)
  extends SamClient(
    "http",
    "bar",
    1,
    checkSubmitWhitelist,
    NoLogging,
    ActorRef.noSender
  )(system, ece, materializer) {

  val authorizedUserCollectionStr: String = "123456789"
  val unauthorizedUserCollectionStr: String = "987654321"

  val notWhitelistedUser: String = "ABC123"

  val userCollectionList: List[Collection] = List(Collection("col1"), Collection("col2"))

  override def collectionsForUser(user: User, httpRequest: HttpRequest): Future[List[Collection]] = {
    val userId = user.userId.value
    if (userId.equalsIgnoreCase(unauthorizedUserCollectionStr)) Future.failed(new Exception(s"Unable to look up collections for user $userId!"))
    else Future.successful(userCollectionList)
  }

  override def requestSubmission(user: User, collection: Collection, cromIamRequest: HttpRequest): Future[Unit] = {
    collection match {
      case c if c.name.equalsIgnoreCase(unauthorizedUserCollectionStr) => Future.failed(SamRegisterCollectionException(StatusCodes.BadRequest))
      case c if c.name.equalsIgnoreCase(authorizedUserCollectionStr) => Future.successful(())
      case _ => Future.successful(())
    }
  }

  override def isSubmitWhitelistedSam(user: User, cromIamRequest: HttpRequest): Future[Boolean] = {
    if (user.userId.value.equalsIgnoreCase(notWhitelistedUser)) Future.successful(false)
    else Future.successful(true)
  }

  override def requestAuth(authorizationRequest: CollectionAuthorizationRequest, cromIamRequest: HttpRequest): Future[Unit] = {
    authorizationRequest.user.userId.value match {
      case `authorizedUserCollectionStr` => Future.successful(())
      case _ => Future.failed(new SamDenialException)
    }
  }
}
