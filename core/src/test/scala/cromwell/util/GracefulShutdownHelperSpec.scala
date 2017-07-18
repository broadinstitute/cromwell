package cromwell.util

import akka.actor.{Actor, ActorLogging, Props}
import akka.testkit.TestProbe
import cats.data.NonEmptyList
import cromwell.core.TestKitSuite
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.scalatest.{FlatSpecLike, Matchers}

class GracefulShutdownHelperSpec extends TestKitSuite with FlatSpecLike with Matchers {
  behavior of "GracefulShutdownHelper"
  
  it should "send ShutdownCommand to actors, wait for them to shutdown, then shut itself down" in {
    val testProbeA = TestProbe()
    val testProbeB = TestProbe()
    
    val testActor = system.actorOf(Props(new Actor with GracefulShutdownHelper with ActorLogging {
      override def receive: Receive = {
        case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.of(testProbeA.ref, testProbeB.ref))
      }
    }))
    
    watch(testActor)

    testActor ! ShutdownCommand

    testProbeA.expectMsg(ShutdownCommand)
    testProbeB.expectMsg(ShutdownCommand)

    // Make sure it's still alive
    expectNoMsg()

    system stop testProbeA.ref

    // Make sure it's still alive
    expectNoMsg()

    system stop testProbeB.ref
    
    expectTerminated(testActor)
  }
}
