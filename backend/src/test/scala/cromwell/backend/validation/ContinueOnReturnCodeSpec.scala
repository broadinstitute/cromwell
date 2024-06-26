package cromwell.backend.validation

import common.assertion.CromwellTimeoutSpec
import org.scalatest.BeforeAndAfterAll
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.wordspec.AnyWordSpecLike

class ContinueOnReturnCodeSpec extends AnyWordSpecLike with CromwellTimeoutSpec with Matchers with BeforeAndAfterAll {
  "Checking for return codes" should {
    "continue on expected return code flags" in {
      val flagTests = Table(("flag", "returnCode", "expectedContinue"),
                            (true, 0, true),
                            (true, 1, true),
                            (false, 0, true),
                            (false, 1, false)
      )

      forAll(flagTests) { (flag, returnCode, expectedContinue) =>
        ContinueOnReturnCodeFlag(flag).continueFor(returnCode) should be(expectedContinue)
      }
    }

    "continue on expected return code sets" in {
      val setTests = Table(
        ("set", "returnCode", "expectedContinue"),
        (Set(0), 0, true),
        (Set(0), 1, false),
        (Set(1), 0, false),
        (Set(1), 1, true),
        (Set(0, 1), 0, true),
        (Set(0, 1), 1, true),
        (Set(0, 1), 2, false)
      )

      forAll(setTests) { (set, returnCode, expectedContinue) =>
        ContinueOnReturnCodeSet(set).continueFor(returnCode) should be(expectedContinue)
      }
    }

    "continue on expected return code string" in {
      val flagTests = Table(("string", "returnCode", "expectedContinue"), ("*", 0, true), ("*", 1, true))

      forAll(flagTests) { (flag, returnCode, expectedContinue) =>
        ContinueOnReturnCodeFlag(flag == "*").continueFor(returnCode) should be(expectedContinue)
      }
    }
  }
}
