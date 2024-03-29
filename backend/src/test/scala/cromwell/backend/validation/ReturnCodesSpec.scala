package cromwell.backend.validation

import common.assertion.CromwellTimeoutSpec
import org.scalatest.BeforeAndAfterAll
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.wordspec.AnyWordSpecLike

class ReturnCodesSpec extends AnyWordSpecLike with CromwellTimeoutSpec with Matchers with BeforeAndAfterAll {
  "Checking for return codes" should {
    "continue on expected return code flags" in {
      val flagTests = Table(("flag", "returnCode", "expectedContinue"), ("*", 0, true), ("*", 1, true))

      forAll(flagTests) { (flag, returnCode, expectedContinue) =>
        ContinueOnReturnCodeFlag(true).continueFor(returnCode) should be(expectedContinue)
      }
    }
  }
}
