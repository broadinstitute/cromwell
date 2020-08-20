package common.util

import cats.effect.IO
import common.util.IORetry.StatefulIoError
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._


class IoRetrySpec extends AnyFlatSpec with Matchers {
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
    val io = IORetry.withRetry(work, 1, Option(3), backoff = Backoff.staticBackoff(10.millis), onRetry = incrementOnRetry)
    val statefulException = the[Exception] thrownBy io.unsafeRunSync()
    statefulException.getCause shouldBe exception
    statefulException.getMessage shouldBe "Attempted 3 times"
  }

}
