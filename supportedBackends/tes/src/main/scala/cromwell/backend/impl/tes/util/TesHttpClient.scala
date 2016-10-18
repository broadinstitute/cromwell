package cromwell.backend.impl.tes.util

import akka.actor.ActorSystem
import spray.http._
import spray.httpx.unmarshalling._
import spray.client.pipelining._
import scala.concurrent.{Future}
import scala.concurrent.ExecutionContext.Implicits.global

object TesClient {

  implicit val system = ActorSystem("tes-http-client")

  def Pipeline[T: FromResponseUnmarshaller]: HttpRequest => Future[T] = sendReceive ~> unmarshal[T]

}
