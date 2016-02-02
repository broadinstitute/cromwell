package cromwell.util

import java.util.UUID

import cromwell.engine.{CromwellFatalException, WorkflowId, WorkflowSourceFiles, WorkflowDescriptor}
import cromwell.logging.WorkflowLogger
import org.scalatest.mock.MockitoSugar
import org.scalatest.{FlatSpec, Matchers}

import scala.language.postfixOps
import scala.concurrent.duration._
import scala.util.Success

class TryUtilSpec extends FlatSpec with Matchers with MockitoSugar {

  class MockWork {
    var counter: Int = _

    def failNTimes(n: Int): Option[Int] => Int = {
      counter = n
      def func(prior: Option[Int]): Int = {
        if (counter > 0) {
          counter -= 1
          throw new IllegalArgumentException("Failed")
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
  }
}
