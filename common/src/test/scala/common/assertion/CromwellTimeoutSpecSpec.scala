package common.assertion

import java.time.OffsetDateTime

import org.scalatest.{Failed, Outcome, Succeeded}
import org.scalatest.exceptions.TestFailedDueToTimeoutException
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.time.Span
import org.scalatest.time.SpanSugar._

class CromwellTimeoutSpecSpec extends AnyFlatSpec with Matchers with CromwellTimeoutSpec {

  override val timeLimit: Span = 5.millis

  override def withFixture(test: NoArgTest): Outcome = {
    val start = OffsetDateTime.now()
    super.withFixture(test) match {
      case Succeeded => Failed("Did not time out")
      case Failed(_: TestFailedDueToTimeoutException) =>
        val timeoutTime = OffsetDateTime.now()
        // Allow 1s of leeway on top of the 5ms timeLimit, but fail if 10 seconds like the Thread.sleep wants:
        val allowableTestSpan = start.plusSeconds(1)
        if (timeoutTime.isBefore(allowableTestSpan)) Succeeded
        else Failed("CromwellTimeoutSpec wrapper did not timeout the thread within 1s (expected timeout after 5ms)")
      case other => Failed(s"Unexpected outcome: $other")
    }
  }

  it should "timeout before reaching an alternative failure state" in {
    Thread.sleep(10000)
    fail("Should have timed out already with a different Exception")
  }
}
