package centaur.api

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class DaemonizedDefaultThreadFactorySpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "DaemonizedDefaultThreadFactory"

  it should "create a non-blocking execution context" in {
    val thread = DaemonizedDefaultThreadFactory.newThread(() => {})
    thread.getName should startWith("daemonpool-thread-")
    thread.isDaemon should be(true)
  }

}
