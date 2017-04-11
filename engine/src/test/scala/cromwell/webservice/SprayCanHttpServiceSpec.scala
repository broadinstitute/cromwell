package cromwell.webservice

import akka.actor.{Actor, Props}
import akka.io.Tcp.Unbound
import akka.testkit._
import akka.util.Timeout
import cromwell.webservice.SprayCanHttpService._
import lenthall.test.actor.TestActorSystem._
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{Assertions, FlatSpec, Matchers}
import spray.can.Http
import spray.routing.HttpService

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Promise}

class SprayCanHttpServiceSpec extends FlatSpec with Matchers with ScalaFutures with Assertions {

  implicit val defaultPatience = PatienceConfig(timeout = Span(5, Seconds), interval = Span(100, Millis))

  import ExecutionContext.Implicits.global

  behavior of "SprayCanHttpService"

  it should "bind and unbind to a random port" in {
    withActorSystem { implicit system =>
      val service = system.actorOf(Props[TestHttpServiceActor])
      implicit val timeout = Timeout(2.seconds.dilated)

      val (port, unbound) = (for {
        boundListener <- service.bind(LoopbackAddress, 0)
        port = boundListener.port
        unbound <- boundListener.unbind()
      } yield (port, unbound)).futureValue

      port should be > 0
      unbound should be(Http.Unbound)

      assert(!system.whenTerminated.isCompleted)
      ()
    }
  }

  it should "bind and unbind to a random port without shutting down" in {
    withActorSystem { implicit system =>
      val service = system.actorOf(Props[TestHttpServiceActor])
      implicit val timeout = Timeout(2.seconds.dilated)

      val (port, unbound) = (for {
        boundListener <- service.bindOrShutdown(LoopbackAddress, 0)
        port = boundListener.port
        unbound <- boundListener.unbind()
      } yield (port, unbound)).futureValue

      port should be > 0
      unbound should be(Http.Unbound)

      assert(!system.whenTerminated.isCompleted)
      ()
    }
  }

  it should "bind and unbind to a random port and then shut down" in {
    withActorSystem { implicit system =>
      val service = system.actorOf(Props[TestHttpServiceActor])
      implicit val timeout = Timeout(2.seconds.dilated)

      val (port, unbound) = (for {
        boundListener <- service.bind(LoopbackAddress, 0)
        port = boundListener.port
        unbound <- boundListener.unbindAndShutdown()
      } yield (port, unbound)).futureValue

      port should be > 0
      unbound should be(Http.Unbound)

      assert(system.whenTerminated.isCompleted)
      ()
    }
  }

  it should "fail to bind to an invalid endpoint" in {
    withActorSystem { implicit system =>
      val service = system.actorOf(Props[TestHttpServiceActor])
      implicit val timeout = Timeout(2.seconds.dilated)

      val exception = (for {
        boundListener <- service.bind(LoopbackAddress, -1)
        _ <- boundListener.unbind()
      } yield Unit).failed.futureValue

      exception should be(an[IllegalArgumentException])
      exception.getMessage should be("port out of range:-1")

      assert(!system.whenTerminated.isCompleted)
      ()
    }
  }

  it should "fail to bind to an invalid endpoint and shutdown" in {
    withActorSystem { implicit system =>
      val service = system.actorOf(Props[TestHttpServiceActor])
      implicit val timeout = Timeout(2.seconds.dilated)

      val exception = (for {
        boundListener <- service.bindOrShutdown(LoopbackAddress, -1)
        _ <- boundListener.unbind()
      } yield Unit).failed.futureValue

      exception should be(an[IllegalArgumentException])
      exception.getMessage should be("port out of range:-1")

      assert(system.whenTerminated.isCompleted)
      ()
    }
  }

  it should "fail to bind to an open endpoint" in {
    withActorSystem { implicit system =>
      val service = system.actorOf(Props[TestHttpServiceActor])
      implicit val timeout = Timeout(2.seconds.dilated)

      val (port, exception, unbound) = (for {
        boundListener <- service.bind(LoopbackAddress, 0)
        port = boundListener.port
        exception <- service.bind(LoopbackAddress, port).failed
        unbound <- boundListener.unbind()
      } yield (port, exception, unbound)).futureValue

      port should be > 0
      exception should be(a[BindFailedException])
      unbound should be(Http.Unbound)

      assert(!system.whenTerminated.isCompleted)
      ()
    }
  }

  it should "return with a bind timeout exception" in {
    withActorSystem { implicit system =>
      val promise = Promise[BoundListener]()
      system.actorOf(Props(classOf[TestSprayCanBindActor], promise, Timeout(100.milliseconds.dilated)))
      promise.future.failed.futureValue should be(a[BindTimeoutException])
      ()
    }
  }

  it should "return with an unexpected bind message exception" in {
    withActorSystem { implicit system =>
      val promise = Promise[BoundListener]()
      val actorRef = system.actorOf(Props(classOf[TestSprayCanBindActor], promise, Timeout(2.seconds.dilated)))
      actorRef ! "unexpected"
      promise.future.failed.futureValue should be(an[UnexpectedBindMessageException])
      ()
    }
  }

  it should "timeout when unbind called twice" in {
    withActorSystem { implicit system =>
      val service = system.actorOf(Props[TestHttpServiceActor])
      implicit val timeout = Timeout(2.seconds.dilated)

      val (port, unbound, exception) = (for {
        boundListener <- service.bind(LoopbackAddress, 0)
        port = boundListener.port
        unbound <- boundListener.unbind()
        exception <- boundListener.unbind().failed // sends a dead letter
      } yield (port, unbound, exception)).futureValue

      port should be > 0
      unbound should be(Http.Unbound)
      exception should be(an[UnbindTimeoutException])

      assert(!system.whenTerminated.isCompleted)
      ()
    }
  }

  it should "timeout when unbind called twice and shutdown" in {
    withActorSystem { implicit system =>
      val service = system.actorOf(Props[TestHttpServiceActor])
      implicit val timeout = Timeout(2.seconds.dilated)

      val (port, unbound, exception) = (for {
        boundListener <- service.bind(LoopbackAddress, 0)
        port = boundListener.port
        unbound <- boundListener.unbind()
        exception <- boundListener.unbindAndShutdown().failed // sends a dead letter
      } yield (port, unbound, exception)).futureValue

      port should be > 0
      unbound should be(Http.Unbound)
      exception should be(an[UnbindTimeoutException])

      assert(system.whenTerminated.isCompleted)
      ()
    }
  }

  it should "return with an unbind failed exception" in {
    withActorSystem { implicit system =>
      val promise = Promise[Unbound]()
      val actorRef = system.actorOf(Props(classOf[TestSprayCanUnbindActor], promise, Timeout(2.seconds.dilated)))
      actorRef ! Http.CommandFailed(Http.Unbind)
      promise.future.failed.futureValue should be(an[UnbindFailedException])
      ()
    }
  }

  it should "return with an unbind timeout exception" in {
    withActorSystem { implicit system =>
      val promise = Promise[Unbound]()
      system.actorOf(Props(classOf[TestSprayCanUnbindActor], promise, Timeout(100.milliseconds.dilated)))
      promise.future.failed.futureValue should be(an[UnbindTimeoutException])
      ()
    }
  }

  it should "return with an unexpected unbind message exception" in {
    withActorSystem { implicit system =>
      val promise = Promise[Unbound]()
      val actorRef = system.actorOf(Props(classOf[TestSprayCanUnbindActor], promise, Timeout(2.seconds.dilated)))
      actorRef ! "unexpected"
      promise.future.failed.futureValue should be(an[UnexpectedUnbindMessageException])
      ()
    }
  }

}

class TestHttpServiceActor extends Actor with HttpService {
  override implicit def actorRefFactory = context

  override def receive = {
    case _ =>
  }
}

class TestSprayCanBindActor(promise: Promise[BoundListener], timeout: Timeout)
  extends SprayCanBindActor(Http.Bind(null, LoopbackAddress), promise, timeout) {
  override def preStart(): Unit = {
    context setReceiveTimeout timeout.duration
  }
}

class TestSprayCanUnbindActor(promise: Promise[Unbound], timeout: Timeout)
  extends SprayCanUnbindActor(null, Http.Unbind(Duration.Zero), promise, timeout) {
  override def preStart(): Unit = {
    context setReceiveTimeout timeout.duration
  }
}
