package cromwell.backend.impl.tes.util

import akka.actor.ActorSystem
import spray.http.HttpRequest
import spray.httpx.unmarshalling._
import spray.client.pipelining._
import scala.concurrent.{Future}

object TesClient {

  implicit val system = ActorSystem("tes-http-client")

  def Pipeline[T: FromResponseUnmarshaller]: HttpRequest => Future[T] = sendReceive ~> unmarshal[T]

}
