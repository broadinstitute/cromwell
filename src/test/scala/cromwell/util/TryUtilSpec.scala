package cromwell.util

import org.scalatest.{FlatSpec, Matchers}

import scala.language.postfixOps
import scala.concurrent.duration._
import scala.util.Success

class TryUtilSpec extends FlatSpec with Matchers {
  def failNTimes(n: Int): Option[Int] => Int = {
    var counter = n
    def func(prior: Option[Int]): Int = {
      if (counter > 0) {
        counter -= 1
        throw new Exception("Failed")
      }
      9
    }
    func
  }

  it should "Retry a function until it works" in {
    val value = TryUtil.retryBlock(
      fn = failNTimes(4),
      retries = Some(5),
      pollingInterval = 50 milliseconds,
      pollingBackOffFactor = 1,
      maxPollingInterval = 10 seconds,
      failMessage = Some(s"failed attempt (on purpose)")
    )
    value shouldEqual Success(9)
  }

  it should "Fail if it hits the max retry count" in {
    val value = TryUtil.retryBlock(
      fn = failNTimes(4),
      retries = Some(4),
      pollingInterval = 50 milliseconds,
      pollingBackOffFactor = 1,
      maxPollingInterval = 10 seconds,
      failMessage = Some(s"failed attempt (on purpose)")
    )
    value.isFailure shouldEqual true
  }
}
