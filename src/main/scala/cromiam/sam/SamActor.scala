package cromiam.sam

import akka.actor.{Actor, ActorLogging, Props}
import cromiam.sam.SamActor._

import scala.concurrent.{ExecutionContext, Promise}
import scala.util.{Failure, Random, Success}
import scala.concurrent.duration._

final class SamActor(userIdHeader: String, allowedUsers: List[String]) extends Actor with ActorLogging {

  implicit val ec: ExecutionContext = context.dispatcher

  override def receive: Receive = {
    case r: SamActorRequest => handleRequest(r)
  }

  private def handleRequest: Function[SamActorRequest, Unit] = {
    case RequestAuthorization(ar, onComplete) =>
      // TODO: Replace with a real IAM lookup:
      val delay = (Random.nextInt(500) + 100).milliseconds
      val message = if (allowedUsers.contains(ar.userIdToken)) Success(()) else Failure(SamDenialException)
      context.system.scheduler.scheduleOnce(delay) { onComplete.complete(message) }
      ()
    case RegisterWorkflow(u, w, onComplete) =>
      // TODO: Replace with a real IAM registration:
      val delay = (Random.nextInt(500) + 100).milliseconds
      context.system.scheduler.scheduleOnce(delay) {
        log.info(s"New workflow ID $w successfully registered to user '$u'")
        onComplete.success(()) }
      ()
  }
}

object SamActor {
  def props(userIdHeader: String, allowedUsers: List[String]): Props = Props(new SamActor(userIdHeader, allowedUsers))

  sealed trait SamActorRequest
  final case class RequestAuthorization(authorizationRequest: AuthorizationRequest, onComplete: Promise[Unit]) extends SamActorRequest
  final case class RegisterWorkflow(user: String, workflowId: String, onComplete: Promise[Unit]) extends SamActorRequest

  case object SamDenialException extends Exception("Access Denied")

  sealed trait AuthorizationRequest {
    def userIdToken: String
    def action: String
  }
  final case class WorkflowAuthorizationRequest(userIdToken: String, workflowId: String, action: String) extends AuthorizationRequest
}
