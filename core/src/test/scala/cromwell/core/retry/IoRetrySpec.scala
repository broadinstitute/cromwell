package cromwell.core.retry

import cats.effect.IO
import cromwell.core.TestKitSuite
import cromwell.core.retry.IORetry.StatefulIoError
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
class IoRetrySpec extends TestKitSuite("retry-spec") with FlatSpecLike with Matchers {
  implicit val timer = IO.timer(ExecutionContext.global)
  implicit val ioError = new StatefulIoError[Int] {
    override def toThrowable(state: Int, throwable: Throwable) = new Exception(s"Attempted $state times", throwable)
  }

  it should "pass state on each retry" in {
    val exception = new Exception("nope")
    var count: Int = 0
    val work = IO {
      if (count == 3) "success"
      else {
        count += 1
        throw exception
      }
    }

    val incrementOnRetry: (Throwable, Int) => Int = (_, s) => s + 1
    val io = IORetry.withRetry(work, 1, Option(3), backoff = SimpleExponentialBackoff(5.millis, 10.millis, 1D), onRetry = incrementOnRetry)
    val statefulException = the[Exception] thrownBy io.unsafeRunSync()
    statefulException.getCause shouldBe exception
    statefulException.getMessage shouldBe "Attempted 3 times"
  }

}
