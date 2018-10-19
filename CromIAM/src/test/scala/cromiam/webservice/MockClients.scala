package cromiam.webservice

import akka.actor.{ActorRef, ActorSystem}
import akka.event.NoLogging
import akka.http.scaladsl.model.{HttpRequest, HttpResponse}
import akka.http.scaladsl.model.StatusCodes._
import akka.stream.ActorMaterializer
import cromiam.auth.{Collection, User}
import cromiam.cromwell.CromwellClient
import cromiam.sam.SamClient

import scala.concurrent.{ExecutionContextExecutor, Future}

final class MockCromwellClient()(implicit system: ActorSystem,
                                 ece: ExecutionContextExecutor,
                                 materializer: ActorMaterializer)
  extends CromwellClient("http", "bar", 1, NoLogging)(system, ece, materializer) {
  override def forwardToCromwell(httpRequest: HttpRequest): Future[HttpResponse] = {
    Future.successful(HttpResponse(status = InternalServerError))
  }
}

class MockSamClient()(implicit system: ActorSystem,
                      ece: ExecutionContextExecutor,
                      materializer: ActorMaterializer)
  extends SamClient("http", "bar", 1, NoLogging, ActorRef.noSender)(system, ece, materializer) {
  override def collectionsForUser(user: User): Future[List[Collection]] = {
    Future.successful(List(Collection("col1"), Collection("col2")))
  }
}
