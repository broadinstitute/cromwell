package cromwell.util

import java.util.UUID

import cromwell.engine.{WorkflowId, WorkflowSourceFiles, WorkflowDescriptor}
import cromwell.logging.WorkflowLogger
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

  val descriptor = WorkflowDescriptor(
    WorkflowId(UUID.randomUUID()),
    WorkflowSourceFiles(
      "workflow w {}",
      "{}",
      "{}"
    )
  )

  val logger = WorkflowLogger("Test", descriptor)

  it should "Retry a function until it works" in {
    val value = TryUtil.retryBlock(
      fn = failNTimes(4),
      retryLimit = Some(5),
      pollingInterval = 50 milliseconds,
      pollingBackOffFactor = 1,
      maxPollingInterval = 10 seconds,
      logger = logger,
      failMessage = Some(s"failed attempt (on purpose)")
    )
    value shouldEqual Success(9)
  }

  it should "Fail if it hits the max retry count" in {
    val value = TryUtil.retryBlock(
      fn = failNTimes(4),
      retryLimit = Some(4),
      pollingInterval = 50 milliseconds,
      pollingBackOffFactor = 1,
      maxPollingInterval = 10 seconds,
      logger = logger,
      failMessage = Some(s"failed attempt (on purpose)")
    )
    value.isFailure shouldEqual true
  }
}
