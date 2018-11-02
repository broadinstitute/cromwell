package cromiam.webservice

import akka.actor.{ActorRef, ActorSystem}
import akka.event.NoLogging
import akka.http.scaladsl.model.{HttpRequest, HttpResponse, StatusCodes}
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
  override def forwardToCromwell(httpRequest: HttpRequest): Future[HttpResponse] = {
    Future.successful(HttpResponse(status = InternalServerError))
  }
}

class MockSamClient()(implicit system: ActorSystem,
                      ece: ExecutionContextExecutor,
                      materializer: ActorMaterializer)
  extends SamClient("http", "bar", 1, NoLogging, ActorRef.noSender)(system, ece, materializer) {

  val unauthorizedUserCollectionStr: String = "123456789"
  val authorizedUserCollectionStr: String = "987654321"

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
