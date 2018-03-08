package cromwell.services.loadcontroller.impl

import akka.actor.ActorRef
import akka.testkit.{TestActorRef, TestProbe}
import cromwell.core.TestKitSuite
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import cromwell.services.loadcontroller.LoadControllerService.{HighLoad, LoadMetric, NormalLoad}
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._

class MemoryLoadControllerSpec extends TestKitSuite with FlatSpecLike with Matchers {
  behavior of "MemoryLoadController"
  
  it should "send memory load metrics" in {
    val registry = TestProbe()
    TestActorRef(new MemoryLoadControllerActorTest(5, 1.second, registry.ref))
    registry.expectMsgPF(2.seconds){
      case metric: LoadMetric => 
        metric.name.head shouldBe "Memory"
        metric.loadLevel shouldBe NormalLoad
    }
  }

  it should "send high memory load metrics iff all load recordings are high" in {
    val registry = TestProbe()
    TestActorRef(new MemoryLoadControllerActorTest(20, 200.milliseconds, registry.ref))
    var metrics: List[LoadMetric] = List.empty
    registry.receiveWhile(3.seconds) {
      case metric: LoadMetric => metrics = metrics :+ metric 
      case _: InstrumentationServiceMessage =>
    }

    // The first 9 should be normal
    metrics.take(9).forall(_.loadLevel == NormalLoad) shouldBe true
    // The rest should be high
    metrics.drop(9).forall(_.loadLevel == HighLoad) shouldBe true 
  }
  
  class MemoryLoadControllerActorTest(threshold: Int, frequency: FiniteDuration, registry: ActorRef) extends MemoryLoadControllerActor(registry) {
    override def getFreeMemory = 10L
    override private[impl] val memoryThreshold = threshold
    override private[impl] val monitoringFrequency = frequency
  }
}
