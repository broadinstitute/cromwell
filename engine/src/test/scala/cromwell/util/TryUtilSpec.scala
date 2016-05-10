package cromwell.util

import cromwell.backend.SimpleExponentialBackoff
import cromwell.core.CromwellFatalException
import cromwell.logging.WorkflowLogger
import org.scalatest.mock.MockitoSugar
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.Success

class TryUtilSpec extends FlatSpec with Matchers with MockitoSugar {

  class TransientException extends Exception
  class MockWork {
    var counter: Int = _

    /**
      * @param n number of times an exception is thrown before succeeding
      * @param transients how many transient exceptions to raise (must be <= n)
      */
    def failNTimes(n: Int, transients: Int = 0): Option[Int] => Int = {
      counter = n
      def func(prior: Option[Int]): Int = {
        if (counter > 0) {
          counter -= 1
          if (counter <= transients) throw new TransientException
          else throw new IllegalArgumentException("Failed")
        }
        9
      }
      func
    }
  }

  val logger = mock[WorkflowLogger]
  val backoff = SimpleExponentialBackoff(50 milliseconds, 10 seconds, 1D)

  it should "Retry a function until it works" in {
    val work = new MockWork
    val value = TryUtil.retryBlock(
      fn = work.failNTimes(4),
      retryLimit = Some(5),
      backoff = backoff,
      logger = logger,
      failMessage = Some(s"failed attempt (on purpose)")
    )
    value shouldEqual Success(9)
    work.counter shouldBe 0
  }

  it should "Fail if it hits the max retry count" in {
    val work = new MockWork
    val value = TryUtil.retryBlock(
      fn = work.failNTimes(4),
      retryLimit = Some(4),
      backoff = backoff,
      logger = logger,
      failMessage = Some(s"failed attempt (on purpose)")
    )
    value.isFailure shouldEqual true
    value.failed.get.isInstanceOf[CromwellFatalException] shouldBe true
  }

  it should "Fail if it hits a fatal exception" in {
    val work = new MockWork
    val value = TryUtil.retryBlock(
      fn = work.failNTimes(4),
      retryLimit = Some(4),
      backoff = backoff,
      logger = logger,
      failMessage = Some(s"failed attempt (on purpose)"),
      isFatal = (t: Throwable) => t.isInstanceOf[IllegalArgumentException]
    )
    value.isFailure shouldEqual true
    value.failed.get.isInstanceOf[CromwellFatalException] shouldBe true
    work.counter shouldBe 3

    val work2 = new MockWork
    val value2 = TryUtil.retryBlock(
      fn = work2.failNTimes(4, 2),
      backoff = backoff,
      retryLimit = Some(4),
      logger = logger,
      failMessage = Some(s"failed attempt (on purpose)"),
      isFatal = (t: Throwable) => t.isInstanceOf[IllegalArgumentException],
      isTransient = (t: Throwable) => t.isInstanceOf[TransientException]
    )

    value2.isFailure shouldEqual true
    value2.failed.get.isInstanceOf[CromwellFatalException] shouldBe true
    work2.counter shouldBe 3
  }

  it should "Not count transient errors against the max limit" in {
    def isTransient(t: Throwable) = t.isInstanceOf[TransientException]

    val work = new MockWork
    val value = TryUtil.retryBlock(
      fn = work.failNTimes(5, 0),
      backoff = backoff,
      retryLimit = Some(5),
      logger = logger,
      failMessage = Some(s"failed attempt (on purpose)"),
      isTransient = isTransient
    )

    value.isFailure shouldEqual true
    work.counter shouldBe 0

    val work2 = new MockWork
    val value2 = TryUtil.retryBlock(
      fn = work2.failNTimes(5, 1),
      backoff = backoff,
      retryLimit = Some(5),
      logger = logger,
      failMessage = Some(s"failed attempt (on purpose)"),
      isTransient = isTransient
    )

    value2.isFailure shouldEqual false
    value2 shouldEqual Success(9)
    work2.counter shouldBe 0

    val work3 = new MockWork
    val value3 = TryUtil.retryBlock(
      fn = work3.failNTimes(6, 4),
      backoff = backoff,
      retryLimit = Some(5),
      logger = logger,
      failMessage = Some(s"failed attempt (on purpose)"),
      isTransient = isTransient
    )

    value3.isFailure shouldEqual false
    value3 shouldEqual Success(9)
    work3.counter shouldBe 0
  }
}
