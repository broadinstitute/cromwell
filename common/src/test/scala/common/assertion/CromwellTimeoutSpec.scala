package common.assertion

import org.scalatest.TestSuite
import org.scalatest.concurrent.{Signaler, ThreadSignaler, TimeLimitedTests}
import org.scalatest.time.Span
import org.scalatest.time.SpanSugar._

trait CromwellTimeoutSpec extends TimeLimitedTests { _: TestSuite =>
  override val timeLimit: Span = 5.minutes
  override val defaultTestSignaler: Signaler = ThreadSignaler
}
