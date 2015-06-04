package cromwell.util

import akka.actor.ActorRef
import akka.pattern.ask
import akka.util.Timeout

import scala.concurrent.duration._
import scala.concurrent.{Await, Future}
import scala.language.{higherKinds, postfixOps, reflectiveCalls}

object ActorUtil {

  implicit val timeout = Timeout(5 seconds)

  /**
   * Performs the following steps:
   *
   * <ol>
   * <li> Sends the specified message to the implicitly passed `ActorRef` via an `ask`.
   * <li> Collects the `Future[Any]` response.
   * <li> Invokes the downcasting function `downcast` on that `Future[Any]`.
   * <li> Issues a blocking `Await.result` on the `Future`.
   * <li> Invokes a `get` on the return of `Await.result`.
   * </ol>
   *
   * The return of `Await.result` is known to support a `get` method: the context bounds for `M` is the "duck type"
   * (formally, a structural type) having a `get` method returning a `U`.
   *
   */
  def messageWaitAndGet[U, M[U] <: {def get : U}](message: AnyRef, downcast: Future[_] => Future[M[U]])(implicit actorRef: ActorRef): U =
    messageAndWait(message, downcast)(actorRef).get

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
