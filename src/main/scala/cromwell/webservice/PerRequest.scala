package cromwell.webservice

import akka.actor.SupervisorStrategy.Stop
import akka.actor.{OneForOneStrategy, _}
import cromwell.webservice.PerRequest._
import spray.http.StatusCodes._
import spray.http._
import spray.httpx.marshalling._
import spray.routing.RequestContext
import spray.httpx.marshalling.ToResponseMarshaller

import scala.concurrent.duration._

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

/*object ImplicitMarshallers {
  implicit val OutputsMarshaller =
    Marshaller.of[(StatusCode, WorkflowOutputResponse)](MediaTypes.`application/json`) { (value, contentType, ctx) =>
      ctx.marshalTo(HttpEntity(contentType, "something"))
    }
  implicit val StatusMarshaller =
    Marshaller.of[(StatusCode, WorkflowStatusResponse)](MediaTypes.`application/json`) { (value, contentType, ctx) =>
      //value._2.toJson.prettyPrint
      ctx.marshalTo(HttpEntity(contentType, "something2"))
    }
  implicit val SubmitMarshaller =
    Marshaller.of[(StatusCode, WorkflowSubmitResponse)](MediaTypes.`application/json`) { (value, contentType, ctx) =>
      ctx.marshalTo(HttpEntity(contentType, "something2"))
    }
  implicit val NoneMarshaller =
    Marshaller.of[(StatusCode, None.type)](MediaTypes.`application/json`) { (value, contentType, ctx) =>
      ctx.marshalTo(HttpEntity(contentType, ""))
    }
}*/

trait PerRequest extends Actor {
  import context._

  def r: RequestContext
  def target: ActorRef
  def message: AnyRef
  def timeout: Duration

  setReceiveTimeout(timeout)
  target ! message

  def receive = {
    case RequestComplete_(response, marshaller) => complete(response)(marshaller)
    case RequestCompleteWithHeaders_(response, headers, marshaller) => complete(response, headers:_*)(marshaller)
    case ReceiveTimeout => complete(GatewayTimeout)
    case x =>
      system.log.error("Unsupported response message sent to PreRequest actor: " + Option(x).getOrElse("null").toString)
      complete(InternalServerError)
  }

  /**
   * Complete the request sending the given response and status code
   * @param response to send to the caller
   * @param marshaller to use for marshalling the response
   * @tparam T the type of the response
   * @return
   */
  private def complete[T](response: T, headers: HttpHeader*)(implicit marshaller: ToResponseMarshaller[T]) = {
    val additionalHeaders = None
    r.withHttpResponseHeadersMapped(h => h ++ headers ++ additionalHeaders).complete(response)
    stop(self)
  }

  override val supervisorStrategy =
    OneForOneStrategy() {
      case e => {
        system.log.error(e, "error processing request: " + r.request.uri)
        r.complete(InternalServerError, e.getMessage)
        Stop
      }
    }
}

object PerRequest {
  sealed trait PerRequestMessage
  /**
   * Report complete, follows same pattern as spray.routing.RequestContext.complete; examples of how to call
   * that method should apply here too. E.g. even though this method has only one parameter, it can be called
   * with 2 where the first is a StatusCode: RequestComplete(StatusCode.Created, response)
   */
  case class RequestComplete[T](response: T)(implicit val marshaller: ToResponseMarshaller[T]) extends PerRequestMessage

  /**
   * Report complete with response headers. To response with a special status code the first parameter can be a
   * tuple where the first element is StatusCode: RequestCompleteWithHeaders((StatusCode.Created, results), header).
   * Note that this is here so that RequestComplete above can behave like spray.routing.RequestContext.complete.
   */
  case class RequestCompleteWithHeaders[T](response: T, headers: HttpHeader*)(implicit val marshaller: ToResponseMarshaller[T]) extends PerRequestMessage

  /** allows for pattern matching with extraction of marshaller */
  private object RequestComplete_ {
    def unapply[T](requestComplete: RequestComplete[T]) = Some((requestComplete.response, requestComplete.marshaller))
  }

  /** allows for pattern matching with extraction of marshaller */
  private object RequestCompleteWithHeaders_ {
    def unapply[T](requestComplete: RequestCompleteWithHeaders[T]) = Some((requestComplete.response, requestComplete.headers, requestComplete.marshaller))
  }

  case class WithProps(r: RequestContext, props: Props, message: AnyRef, timeout: Duration) extends PerRequest {
    lazy val target = context.actorOf(props)
  }
}

/**
 * Provides factory methods for creating per request actors
 */
trait PerRequestCreator {
  implicit def actorRefFactory: ActorRefFactory

  def perRequest(r: RequestContext, props: Props, message: AnyRef, timeout: Duration = 1 minutes) =
    actorRefFactory.actorOf(Props(new WithProps(r, props, message, timeout)))
}