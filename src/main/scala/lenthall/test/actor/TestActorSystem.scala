package lenthall.test.actor

import akka.actor.ActorSystem

import scala.concurrent.Await
import scala.concurrent.duration.Duration

/**
 * Provides an individual ActorSystem for tests.
 * An alternative to the NOT thread-safe TestKit, that may run into problems with parallel testing.
 * http://doc.akka.io/api/akka/2.3.14/index.html#akka.testkit.TestKit
 */
object TestActorSystem {
  def withActorSystem(f: ActorSystem => Unit): Unit = {
    withActorSystem("test-system", awaitTermination = true)(f)
  }

  def withActorSystem(name: String, awaitTermination: Boolean)(f: ActorSystem => Unit): Unit = {
    val system = ActorSystem(name)
    try {
      f(system)
    } finally {
      system.terminate()
      if (awaitTermination)
        Await.ready(system.whenTerminated, Duration.Inf)
    }
  }
}
