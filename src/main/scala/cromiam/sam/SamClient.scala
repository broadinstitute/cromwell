package cromiam.sam

import akka.actor.{ActorRef, ActorRefFactory}
import cromiam.sam.SamActor.{AuthorizationRequest, RegisterWorkflow, RequestAuthorization}

import scala.concurrent.{ExecutionContext, Future, Promise}
import cats.implicits._

trait SamClient {
  protected def samActor: ActorRef
  protected def actorRefFactory: ActorRefFactory

  protected def requestAuth(authorizationRequest: AuthorizationRequest): Future[Unit] = {
    val completeMe: Promise[Unit] = Promise()
    samActor.tell(RequestAuthorization(authorizationRequest, completeMe), ActorRef.noSender)
    completeMe.future
  }

  protected def registerCreation(user: String, workflowIds: List[String])(implicit ec: ExecutionContext): Future[Unit] = {
    (workflowIds traverse { workflowId =>
      val completeMe: Promise[Unit] = Promise()
      samActor.tell(RegisterWorkflow(user, workflowId, completeMe), ActorRef.noSender)
      completeMe.future
    }).void
  }
}
