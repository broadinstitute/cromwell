package lenthall.test.actor

import akka.actor.ActorSystem

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
      system.shutdown()
      if (awaitTermination)
        system.awaitTermination()
    }
  }
}
