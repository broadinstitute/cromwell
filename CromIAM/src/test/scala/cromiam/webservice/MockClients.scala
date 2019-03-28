package cromiam.webservice

import akka.actor.{ActorRef, ActorSystem}
import akka.event.NoLogging
import akka.http.scaladsl.model.StatusCodes._
import akka.http.scaladsl.model.{HttpRequest, HttpResponse, StatusCodes}
import akka.stream.ActorMaterializer
import cats.effect.IO
import cromiam.auth.{Collection, User}
import cromiam.cromwell.CromwellClient
import cromiam.sam.SamClient
import cromiam.sam.SamClient.{CollectionAuthorizationRequest, SamDenialException, SamRegisterCollectionException}
import cromwell.api.model._

import scala.concurrent.ExecutionContextExecutor

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

  def checkHeaderAuthorization(httpRequest: HttpRequest): FailureResponseOrT[HttpResponse] = {
    val userIdHeader = httpRequest.headers.find(header => header.name.equalsIgnoreCase("OIDC_CLAIM_user_id"))

    userIdHeader match {
      case Some(header) =>
        header.value match {
          case headerValue if headerValue.equalsIgnoreCase(authorizedUserCollectionStr) =>
            FailureResponseOrT.pure(HttpResponse(status = OK, entity = "Response from Cromwell"))
          case headerValue =>
            val message = s"This is unexpected! Cromwell should not receive request from unauthorized user! " +
              s"OIDC_CLAIM_user_id: $headerValue is not authorized."
            FailureResponseOrT[IO, HttpResponse, HttpResponse](IO.raiseError(new RuntimeException(message)))
        }
      case None =>
        val message = "This is unexpected! No OIDC_CLAIM_user_id provided for authorization!"
        FailureResponseOrT[IO, HttpResponse, HttpResponse](IO.raiseError(new RuntimeException(message)))
    }
  }

  override def forwardToCromwell(httpRequest: HttpRequest): FailureResponseOrT[HttpResponse] = {
    val versionRoutePath = s"/engine/$version/version"
    val womtoolRoutePath = s"/api/womtool/$version/describe"

    httpRequest.uri.path.toString match {
      //version endpoint doesn't require authentication
      case `versionRoutePath` =>
        FailureResponseOrT.pure(HttpResponse(status = OK, entity = "Response from Cromwell"))
      // womtool endpoint requires authn which it gets for free from the proxy, does not care about authz
      case `womtoolRoutePath` =>
        FailureResponseOrT.pure(HttpResponse(status = OK, entity = "Hey there, workflow describer"))
      case _ => checkHeaderAuthorization(httpRequest)
    }
  }

  override def getRootWorkflow(workflowId: String,
                               user: User,
                               cromIamRequest: HttpRequest): FailureResponseOrT[String] = {
    workflowId match {
      case `subworkflowId` | `rootWorkflowIdWithCollection` =>
        FailureResponseOrT.pure(rootWorkflowIdWithCollection)
      case `anotherRootWorkflowIdWithCollection` => FailureResponseOrT.pure(anotherRootWorkflowIdWithCollection)
      case _ => FailureResponseOrT.pure(workflowIdWithoutCollection)
    }
  }

  override def collectionForWorkflow(workflowId: String,
                                     user: User,
                                     cromIamRequest: HttpRequest): FailureResponseOrT[Collection] = {
    workflowId match {
      case `rootWorkflowIdWithCollection` | `anotherRootWorkflowIdWithCollection` =>
        FailureResponseOrT.pure(userCollection)
      case _ =>
        val exception = new IllegalArgumentException(s"Workflow $workflowId has no associated collection")
        FailureResponseOrT.left(IO.raiseError[HttpResponse](exception))
    }
  }
}

class MockSamClient(checkSubmitWhitelist: Boolean = true,
                    samResponseOption: Option[HttpResponse] = None)
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

  override def collectionsForUser(user: User, httpRequest: HttpRequest): FailureResponseOrT[List[Collection]] = {
    samResponseOption match {
      case Some(samResponse) => FailureResponseOrT.left(IO.pure(samResponse))
      case None =>
        val userId = user.userId.value
        if (userId.equalsIgnoreCase(unauthorizedUserCollectionStr)) {
          val exception = new Exception(s"Unable to look up collections for user $userId!")
          FailureResponseOrT.left(IO.raiseError[HttpResponse](exception))
        } else {
          FailureResponseOrT.pure(userCollectionList)
        }
    }
  }

  override def requestSubmission(user: User,
                                 collection: Collection,
                                 cromIamRequest: HttpRequest): FailureResponseOrT[Unit] = {
    samResponseOption match {
      case Some(samResponse) => FailureResponseOrT.left(IO.pure(samResponse))
      case None =>
        collection match {
          case c if c.name.equalsIgnoreCase(unauthorizedUserCollectionStr) =>
            val exception = SamRegisterCollectionException(StatusCodes.BadRequest)
            FailureResponseOrT.left(IO.raiseError[HttpResponse](exception))
          case c if c.name.equalsIgnoreCase(authorizedUserCollectionStr) => FailureResponseOrT.pure(())
          case _ => FailureResponseOrT.pure(())
        }
    }
  }

  override def isSubmitWhitelistedSam(user: User, cromIamRequest: HttpRequest): FailureResponseOrT[Boolean] = {
    samResponseOption match {
      case Some(samResponse) => FailureResponseOrT.left(IO.pure(samResponse))
      case None => FailureResponseOrT.pure(!user.userId.value.equalsIgnoreCase(notWhitelistedUser))
    }
  }

  override def requestAuth(authorizationRequest: CollectionAuthorizationRequest,
                           cromIamRequest: HttpRequest): FailureResponseOrT[Unit] = {
    samResponseOption match {
      case Some(samResponse) => FailureResponseOrT.left(IO(samResponse))
      case None =>
        authorizationRequest.user.userId.value match {
          case `authorizedUserCollectionStr` => FailureResponseOrT.pure(())
          case _ => FailureResponseOrT.left(IO.raiseError[HttpResponse](new SamDenialException))
        }
    }
  }
}
