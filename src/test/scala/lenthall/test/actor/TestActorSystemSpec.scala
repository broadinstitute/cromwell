package lenthall.test.actor

import akka.actor.ActorSystem
import akka.testkit._
import lenthall.test.actor.TestActorSystem._
import org.scalatest.{Assertions, FlatSpec, Matchers}

import scala.concurrent.Await
import scala.concurrent.duration._

class TestActorSystemSpec extends FlatSpec with Matchers with Assertions {
  behavior of "TestActorSystem"

  it should "create and shutdown a default named system" in {
    var lentSystem: ActorSystem = null
    withActorSystem { system =>
      system.name should be("test-system")
      lentSystem = system
      assert(!system.whenTerminated.isCompleted)
    }
    // Should always be completed
    assert(lentSystem.whenTerminated.isCompleted)
  }

  it should "create and start shutting down a named system" in {
    var timeout: Duration = null
    var lentSystem: ActorSystem = null
    withActorSystem("my-name", awaitTermination = false) { system =>
      system.name should be("my-name")
      timeout = 1.second.dilated(system)
      lentSystem = system
      assert(!system.whenTerminated.isCompleted)
    }
    // As we told the block not to wait, we're going to wait a bit for this termination
    Await.ready(lentSystem.whenTerminated, timeout)
  }
}
