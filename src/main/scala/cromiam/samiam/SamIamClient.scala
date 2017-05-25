package cromiam.samiam

import akka.actor.{ActorRef, ActorRefFactory}
import cromiam.samiam.SamIamActor.{AuthorizationRequest, RegisterWorkflow, RequestAuthorization}

import scala.concurrent.{ExecutionContext, Future, Promise}
import cats.implicits._

trait SamIamClient {
  protected def samIamActor: ActorRef
  protected def actorRefFactory: ActorRefFactory

  protected def requestAuth(authorizationRequest: AuthorizationRequest): Future[Unit] = {
    val completeMe: Promise[Unit] = Promise()
    samIamActor.tell(RequestAuthorization(authorizationRequest, completeMe), ActorRef.noSender)
    completeMe.future
  }

  protected def registerCreation(user: String, workflowIds: List[String])(implicit ec: ExecutionContext): Future[Unit] = {
    (workflowIds traverse { workflowId =>
      val completeMe: Promise[Unit] = Promise()
      samIamActor.tell(RegisterWorkflow(user, workflowId, completeMe), ActorRef.noSender)
      completeMe.future
    }).void
  }
}
