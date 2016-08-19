package centaur.api

import akka.actor.ActorSystem
import centaur.CentaurConfig
import spray.client.pipelining._
import spray.http.HttpRequest
import spray.httpx.unmarshalling._

import scala.concurrent.{Await, Future}
import scala.concurrent.duration.FiniteDuration
import scala.concurrent.ExecutionContext.Implicits.global
import scala.util.Try

object CromwellClient {
  // Spray needs an implicit ActorSystem
  implicit val system = ActorSystem("centaur-foo")
  def Pipeline[T: FromResponseUnmarshaller]: HttpRequest => Future[T] = sendReceive ~> unmarshal[T]

  /**
    * Ensure that the Future completes within the specified timeout. If it does not, or if the Future fails,
    * will return a Failure, otherwise a Success
    */
  def awaitFutureCompletion[T](x: Future[T], timeout: FiniteDuration) = Try(Await.result(x, timeout))
  def sendReceiveFutureCompletion[T](x: Future[T]) = awaitFutureCompletion(x, CentaurConfig.sendReceiveTimeout)
}
