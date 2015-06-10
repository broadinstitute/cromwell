package cromwell.util

import akka.actor.ActorRef
import akka.pattern.ask
import akka.util.Timeout

import scala.concurrent.{Await, Future}
import scala.concurrent.duration._
import scala.language.postfixOps

object ActorTestUtil {

  implicit val timeout = Timeout(5 seconds)

  /**
   * Performs the following steps:
   *
   * <ol>
   * <li> Sends the specified message to the implicitly passed `ActorRef` via an `ask`.
   * <li> Collects the `Future[Any]` response.
   * <li> Invokes the downcasting function `downcast` on that `Future[Any]`.
   * <li> Issues a blocking `Await.result` on the `Future`.
   * </ol>
   *
   */
  def messageAndWait[M](message: AnyRef, downcast: Future[_] => Future[M])(implicit actorRef: ActorRef): M = {
    val futureAny = actorRef ? message
    Await.result(downcast(futureAny), 5 seconds)
  }
}
