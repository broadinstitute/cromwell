package lenthall.spray

import java.net.{InetAddress, InetSocketAddress}

import akka.actor._
import akka.io.Tcp.Unbound
import akka.io.{IO, Inet}
import akka.util.Timeout
import spray.can.Http
import spray.can.Http.Bound
import spray.can.server.ServerSettings
import spray.io.ServerSSLEngineProvider

import scala.collection.immutable
import scala.concurrent._
import scala.concurrent.duration.Duration

/**
 * Adds bind and unbind to http services running on spray-can.
 *
 * {{{
 * val httpServiceActorRef: ActorRef = ...
 * httpServiceActorRef.bind(LoopbackAddress, 8080) onSuccess {
 *   case boundListener =>
 *     boundListener.unbind() // Unbind the listener
 * }
 * }}}
 */
object SprayCanHttpService {

  val LoopbackAddress = InetAddress.getLoopbackAddress.getHostAddress
  val AnyAddress = new InetSocketAddress(0).getAddress.getHostAddress

  implicit class EnhancedSprayCanActor(val service: ActorRef) extends AnyVal {

    def bind(interface: String = LoopbackAddress, port: Int = 80, backlog: Int = 100,
             options: immutable.Traversable[Inet.SocketOption] = Nil, settings: Option[ServerSettings] = None)
            (implicit ec: ExecutionContext, timeout: Timeout, actorSystem: ActorSystem,
             sslEngineProvider: ServerSSLEngineProvider): Future[BoundListener] = {
      Future(Http.Bind(service, interface, port, backlog, options, settings)) flatMap { httpBind =>
        val promise = Promise[BoundListener]()
        actorSystem.actorOf(Props(classOf[SprayCanBindActor], httpBind, promise, timeout))
        promise.future
      }
    }

    def bindOrShutdown(interface: String = LoopbackAddress, port: Int = 80, backlog: Int = 100,
                       options: immutable.Traversable[Inet.SocketOption] = Nil, settings: Option[ServerSettings] = None)
                      (implicit ec: ExecutionContext, actorSystem: ActorSystem, timeout: Timeout,
                       sslEngineProvider: ServerSSLEngineProvider): Future[BoundListener] = {
      bind(interface, port, backlog, options, settings) recover {
        case throwable =>
          actorSystem.log.error(throwable, s"Binding failed interface $interface port $port")
          actorSystem.terminate()
          /*
          Using scala's Await.ready due to akka recommendations:
          - http://doc.akka.io/docs/akka/2.4/project/migration-guide-2.3.x-2.4.x.html#Actor_system_shutdown
          - https://github.com/akka/akka/blob/v2.4.10/akka-actor/src/main/scala/akka/actor/ActorSystem.scala#L664-L665
           */
          Await.ready(actorSystem.whenTerminated, Duration.Inf)
          throw throwable
      }
    }

  }

  case class BoundListener(bound: Bound, listener: ActorRef) {

    def unbind(duration: Duration = Duration.Zero)
              (implicit ec: ExecutionContext, actorSystem: ActorSystem, timeout: Timeout): Future[Unbound] = {
      Future(Http.Unbind(duration)) flatMap { httpUnbind =>
          val promise = Promise[Unbound]()
          actorSystem.actorOf(Props(classOf[SprayCanUnbindActor], listener, httpUnbind, promise, timeout))
          promise.future
      }
    }

    def unbindAndShutdown(duration: Duration = Duration.Zero)
                         (implicit ec: ExecutionContext, actorSystem: ActorSystem, timeout: Timeout)
    : Future[Unbound] = {
      unbind(duration) andThen {
        case _ =>
          actorSystem.terminate()
          /*
          Using scala's Await.ready due to akka recommendations:
          - http://doc.akka.io/docs/akka/2.4/project/migration-guide-2.3.x-2.4.x.html#Actor_system_shutdown
          - https://github.com/akka/akka/blob/v2.4.10/akka-actor/src/main/scala/akka/actor/ActorSystem.scala#L664-L665
           */
          Await.ready(actorSystem.whenTerminated, Duration.Inf)
      }
    }

    // Alias for testing
    private[spray] def port = bound.localAddress.getPort

  }

  trait BindException {
    this: Exception =>
    val httpBind: Http.Bind
  }

  class BindFailedException(message: String, val httpBind: Http.Bind)
    extends RuntimeException(message) with BindException

  class BindTimeoutException(message: String, val httpBind: Http.Bind, val timeout: Timeout)
    extends TimeoutException(message) with BindException

  class UnexpectedBindMessageException(message: String, val httpBind: Http.Bind, val unexpected: Any)
    extends RuntimeException(message) with BindException

  class SprayCanBindActor(httpBind: Http.Bind, promise: Promise[BoundListener], timeout: Timeout) extends Actor {

    override def preStart() {
      IO(Http)(context.system) ! httpBind
      context setReceiveTimeout timeout.duration
    }

    override def receive = {
      case bound: Bound =>
        val sentBy = sender()
        promise trySuccess BoundListener(bound, sentBy)
        context stop self
      case failed: Http.CommandFailed =>
        promise tryFailure new BindFailedException(s"Failed to bind to ${httpBind.endpoint}: $failed", httpBind)
        context stop self
      case ReceiveTimeout =>
        promise tryFailure new BindTimeoutException(
          s"Timeout ${timeout.duration} during bind to ${httpBind.endpoint}", httpBind, timeout)
        context stop self
      case unexpected =>
        promise tryFailure new UnexpectedBindMessageException(
          s"Unexpected message during bind to ${httpBind.endpoint}: $unexpected", httpBind, unexpected)
        context stop self
    }
  }

  trait UnbindException {
    this: Exception =>
    val httpUnbind: Http.Unbind
  }

  class UnbindFailedException(message: String, val httpUnbind: Http.Unbind)
    extends RuntimeException(message) with UnbindException

  class UnbindTimeoutException(message: String, val httpUnbind: Http.Unbind, val timeout: Timeout)
    extends TimeoutException(message) with UnbindException

  class UnexpectedUnbindMessageException(message: String, val httpUnbind: Http.Unbind, val unexpected: Any)
    extends RuntimeException(message) with UnbindException

  class SprayCanUnbindActor(listener: ActorRef, httpUnbind: Http.Unbind, promise: Promise[Unbound], timeout: Timeout)
    extends Actor {

    override def preStart() {
      listener ! httpUnbind
      context setReceiveTimeout timeout.duration
    }

    override def receive = {
      case unbound: Unbound =>
        promise trySuccess unbound
        context stop self
      case failed: Http.CommandFailed =>
        promise tryFailure new UnbindFailedException(s"Failed to unbind: $failed", httpUnbind)
        context stop self
      case ReceiveTimeout =>
        promise tryFailure new UnbindTimeoutException(s"Timeout ${timeout.duration} during unbind", httpUnbind, timeout)
        context stop self
      case unexpected =>
        promise tryFailure new UnexpectedUnbindMessageException(
          s"Unexpected message during unbind: $unexpected", httpUnbind, unexpected)
        context stop self
    }
  }

}
