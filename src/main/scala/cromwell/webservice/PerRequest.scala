package cromwell.webservice

import akka.actor._
import akka.actor.SupervisorStrategy.Stop
import cromwell.webservice.PerRequest._
import spray.http.StatusCodes._
import spray.routing.RequestContext
import akka.actor.OneForOneStrategy
import spray.httpx.Json4sSupport
import scala.concurrent.duration._
import org.json4s.DefaultFormats
import spray.http.StatusCode

/**
 * This actor controls the lifecycle of a request. It is responsible for forwarding the initial message
 * to a target handling actor. This actor waits for the target actor signal completion (via a message),
 * timeout, or handle an exception. It is this actors responsibility to respond to the request and
 * shutdown itself and child actors.
 *
 * Request completion can be signaled in 2 ways:
 * 1) with just a response object
 * 2) with a RequestComplete message which can specify http status code as well as the response
 */
trait PerRequest extends Actor with Json4sSupport {

  import context._

  val json4sFormats = DefaultFormats

  def r: RequestContext
  def target: ActorRef
  def message: AnyRef
  def timeout: Duration

  setReceiveTimeout(timeout)
  target ! message

  def receive = {
    case RequestComplete(statusCode, response) => complete(statusCode, response)
    case ReceiveTimeout => complete(GatewayTimeout, "Request timeout")
    case response: AnyRef => complete(OK, response)
  }

  def complete[T <: AnyRef](status: StatusCode, obj: T) = {
    r.complete(status, obj)
    stop(self)
  }

  override val supervisorStrategy =
    OneForOneStrategy() {
      case e => {
        complete(InternalServerError, e.getMessage)
        Stop
      }
    }
}

object PerRequest {
  /**
   * Message used to control the http status code sent to client
   */
  case class RequestComplete(statusCode: StatusCode, response: AnyRef)

  case class WithProps(r: RequestContext, props: Props, message: AnyRef, timeout: Duration) extends PerRequest {
    lazy val target = context.actorOf(props)
  }
}

/**
 * Provides factory methods for creating per request actors
 */
trait PerRequestCreator {
  implicit def actorRefFactory: ActorRefFactory

  def perRequest(r: RequestContext, props: Props, message: AnyRef, timeout: Duration = 1.minutes) =
    actorRefFactory.actorOf(Props(new WithProps(r, props, message, timeout)))
}