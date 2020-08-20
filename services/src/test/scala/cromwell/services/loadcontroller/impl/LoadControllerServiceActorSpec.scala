package cromwell.services.loadcontroller.impl

import akka.actor.Kill
import akka.routing.Listen
import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import cromwell.core.TestKitSuite
import cromwell.services.loadcontroller.LoadControllerService.{HighLoad, LoadMetric, NormalLoad}
import cromwell.services.loadcontroller.impl.LoadControllerServiceActor.ActorAndMetric
import cromwell.services.loadcontroller.impl.LoadControllerServiceActorSpec._
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import scala.concurrent.duration._

object LoadControllerServiceActorSpec {
  val Config = ConfigFactory.parseString("control-frequency = 1 second")
}

class LoadControllerServiceActorSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with Eventually with ImplicitSender {
  behavior of "LoadControllerServiceActor"
  override implicit val patienceConfig: PatienceConfig = PatienceConfig(scaled(2.seconds))

  it should "record metrics" in {
    val loadActor = TestActorRef(new LoadControllerServiceActor(Config, Config, TestProbe().ref))
    loadActor ! LoadMetric("itsme", NormalLoad)
    loadActor ! LoadMetric("itsotherme", HighLoad)
    loadActor.underlyingActor.loadMetrics(ActorAndMetric(self, NonEmptyList.one("itsme"))) shouldBe NormalLoad
    loadActor.underlyingActor.loadMetrics(ActorAndMetric(self, NonEmptyList.one("itsotherme"))) shouldBe HighLoad
  }

  it should "update global load level periodically" in {
    val loadActor = TestActorRef(new LoadControllerServiceActor(Config, Config, TestProbe().ref))
    loadActor.underlyingActor.loadLevel shouldBe NormalLoad
    loadActor ! LoadMetric("itsme", NormalLoad)
    loadActor ! LoadMetric("itsotherme", HighLoad)
    eventually {
      loadActor.underlyingActor.loadLevel shouldBe HighLoad
    }
    loadActor ! LoadMetric("itsotherme", NormalLoad)
    eventually {
      loadActor.underlyingActor.loadLevel shouldBe NormalLoad
    }
  }

  it should "evict monitored actors from the metric store if they die" in {
    val snd = TestProbe().ref
    val loadActor = TestActorRef(new LoadControllerServiceActor(Config, Config, TestProbe().ref))
    // Send 2 metrics from the same actor
    loadActor.tell(LoadMetric("itsme", NormalLoad), snd)
    loadActor.tell(LoadMetric("itsmeagain", HighLoad), snd)

    // Check that we've got them
    loadActor.underlyingActor.loadMetrics(ActorAndMetric(snd, NonEmptyList.one("itsme"))) shouldBe NormalLoad
    loadActor.underlyingActor.loadMetrics(ActorAndMetric(snd, NonEmptyList.one("itsmeagain"))) shouldBe HighLoad

    // Check that we're monitoring this actor
    loadActor.underlyingActor.monitoredActors shouldBe Set(snd)

    // Kill the sender
    snd ! Kill

    // Check that it's gone
    eventually {
      loadActor.underlyingActor.loadMetrics shouldBe empty
      loadActor.underlyingActor.monitoredActors shouldBe empty
    }
  }

  it should "send updates to listeners when necessary" in {
    val loadActor = TestActorRef(new LoadControllerServiceActor(Config, Config, TestProbe().ref))
    loadActor ! Listen(self)

    // Raise the load level
    loadActor ! LoadMetric("itsme", HighLoad)

    // We should get an alert
    expectMsg(max = 2.seconds, HighLoad)

    // And only one
    expectNoMessage(2.seconds)

    // Set the load level back to normal
    loadActor ! LoadMetric("itsme", NormalLoad)

    // We should get an alert that the load is back to normal
    expectMsg(max = 2.seconds, NormalLoad)
  }
}
