package cromwell.backend.impl.tes.util

import akka.actor.ActorSystem
import spray.client.pipelining._
import spray.http.HttpRequest
import spray.httpx.unmarshalling._

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future

object TesHttpClient {

  implicit val system = ActorSystem("tes-http-client")

  def Pipeline[T: FromResponseUnmarshaller]: HttpRequest => Future[T] = sendReceive ~> unmarshal[T]

}
