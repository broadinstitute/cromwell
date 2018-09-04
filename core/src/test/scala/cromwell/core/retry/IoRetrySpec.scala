package cromwell.core.retry

import cats.Show
import cats.effect.IO
import cromwell.core.TestKitSuite
import cromwell.core.retry.IORetry.StatefulIoException
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
class IoRetrySpec extends TestKitSuite("retry-spec") with FlatSpecLike with Matchers {
  implicit val timer = IO.timer(ExecutionContext.global)
  implicit val showState: Show[Int] = (t: Int) => s"Attempted $t times"

  it should "count attempts when failed" in {
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
    val statefulException = the[StatefulIoException[Int]] thrownBy io.unsafeRunSync()
    statefulException.cause shouldBe exception
    statefulException.state shouldBe 3
  }

}
