package cromiam.webservice

import akka.actor.{ActorRef, ActorSystem}
import akka.event.NoLogging
import akka.http.scaladsl.model.{HttpHeader, HttpRequest, HttpResponse, StatusCodes}
import akka.http.scaladsl.model.StatusCodes._
import akka.stream.ActorMaterializer
import cromiam.auth.{Collection, User}
import cromiam.cromwell.CromwellClient
import cromiam.sam.SamClient
import cromiam.sam.SamClient.SamRegisterCollectionException

import scala.concurrent.{ExecutionContextExecutor, Future}

final class MockCromwellClient()(implicit system: ActorSystem,
                                 ece: ExecutionContextExecutor,
                                 materializer: ActorMaterializer)
  extends CromwellClient("http", "bar", 1, NoLogging, ActorRef.noSender)(system, ece, materializer) {

  val authorizedUserCollectionStr: String = "123456789"

  def checkHeaderAuthorization(header: HttpHeader): Future[HttpResponse] = {
    header.value match {
      case s if s.equalsIgnoreCase(authorizedUserCollectionStr) => Future.successful(HttpResponse(status = OK))
      case s => Future.failed(new Exception(s"This is unexpected! Cromwell should not receive request from unauthorized user! OIDC_CLAIM_user_id: $s is not authorized."))
    }
  }

  override def forwardToCromwell(httpRequest: HttpRequest): Future[HttpResponse] = {
    val userIdHeader = httpRequest.headers.find(header => header.name.equalsIgnoreCase("OIDC_CLAIM_user_id"))

    userIdHeader match {
      case Some(header) => checkHeaderAuthorization(header)
      case None => Future.failed(new Exception("No OIDC_CLAIM_user_id provided for authorization!"))
    }
  }
}

class MockSamClient()(implicit system: ActorSystem,
                      ece: ExecutionContextExecutor,
                      materializer: ActorMaterializer)
  extends SamClient("http", "bar", 1, NoLogging, ActorRef.noSender)(system, ece, materializer) {

  val authorizedUserCollectionStr: String = "123456789"
  val unauthorizedUserCollectionStr: String = "987654321"

  override def collectionsForUser(user: User, httpRequest: HttpRequest): Future[List[Collection]] = {
    Future.successful(List(Collection("col1"), Collection("col2")))
  }

  override def requestSubmission(user: User, collection: Collection, cromIamRequest: HttpRequest): Future[Unit] = {
    collection match {
      case c if c.name.equalsIgnoreCase(unauthorizedUserCollectionStr) => Future.failed(SamRegisterCollectionException(StatusCodes.BadRequest))
      case c if c.name.equalsIgnoreCase(authorizedUserCollectionStr) => Future.successful(())
      case _ => Future.successful(())
    }
  }

  override def isSubmitWhitelisted(user: User, cromIamRequest: HttpRequest): Future[Boolean] = {
    Future.successful(true)
  }
}
