package cromwell.core.retry

import cromwell.core.retry.Retry._
import cromwell.core.{CromwellFatalException, TestKitSuite}
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.Future

class RetrySpec extends TestKitSuite("retry-spec") with FlatSpecLike with Matchers with ScalaFutures {
  class TransientException extends Exception
  class MockWork(n: Int, transients: Int = 0) {
    implicit val ec = system.dispatcher

    var counter = n

    def doIt(): Future[Int] = {
      if (counter == 0)
        Future.successful(9)
      else {
        counter -= 1
        val ex = if (counter <= transients) new TransientException else new IllegalArgumentException("Failed")
        Future.failed(ex)
      }
    }
  }

  implicit val defaultPatience = PatienceConfig(timeout = Span(30, Seconds), interval = Span(100, Millis))

  private def runRetry(retries: Int,
                       work: MockWork,
                       isTransient: Throwable => Boolean = Retry.throwableToFalse,
                       isFatal: Throwable => Boolean = Retry.throwableToFalse): Future[Int] = {

    withRetry(
      f = work.doIt,
      maxRetries = Option(retries),
      isTransient = isTransient,
      isFatal = isFatal
    )
  }

  "Retry" should "retry a function until it works" in {
    val work = new MockWork(2)

    whenReady(runRetry(3, work)) { x =>
      x shouldBe 9
      work.counter shouldBe 0
    }
  }

  it should "fail if it hits the max retry count" in {
    whenReady(runRetry(1, new MockWork(3)).failed) { x =>
      x shouldBe an [CromwellFatalException]
    }
  }

  it should "fail if it hits a fatal exception" in {
    val work = new MockWork(3)

    whenReady(runRetry(3, work, isFatal = (t: Throwable) => t.isInstanceOf[IllegalArgumentException]).failed) { x =>
      x shouldBe an [CromwellFatalException]
      work.counter shouldBe 2
    }

    val work2 = new MockWork(4, 2)
    val retry = runRetry(4,
      work2,
      isFatal = (t: Throwable) => t.isInstanceOf[IllegalArgumentException],
      isTransient = (t: Throwable) => t.isInstanceOf[TransientException])

    whenReady(retry.failed) { x =>
      x shouldBe an [CromwellFatalException]
      work2.counter shouldBe 3
    }
  }

  it should "not count transient errors against the max limit" in {
    val work = new MockWork(3, 1)
    whenReady(runRetry(3, work, isTransient = (t: Throwable) => t.isInstanceOf[TransientException])) { x =>
      x shouldBe 9
      work.counter shouldBe 0
    }
  }
}
