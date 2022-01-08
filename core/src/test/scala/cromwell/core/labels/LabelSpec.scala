package cromwell.core.labels

import cats.data.Validated.{Invalid, Valid}
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class LabelSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "Labels"

  /**
    * In the format 'to validate', 'expected result'
    */
  val goodLabelKeys = List(
    "cromwell-root-workflow-id",
    "cromwell-11f2468c-39d6-4be3-85c8-32735c01e66b",
    "just-the-right-length-just-the-right-length-just-the-right-leng",
    "11f2468c-39d6-4be3-85c8-32735c01e66b",
    "0-cromwell-root-workflow-id",
    "cromwell-root-workflow-id-",
    "0-cromwell-root-workflow-id-",
    "Cromwell-root-workflow-id",
    "now valid 255 character key-at vero eos et accusamus et iusto odio dignissimos ducimus qui blanditiis praesentium voluptatum deleniti atque corrupti quos dolores et quas molestias excepturi sint occaecati cupiditate non provident, similique sunt in culpa",
    "!@#$%^&*()_+={}[]:;'<>?,./`~"
  )

  val goodLabelValues = List(
    "11f2468c-39d6-4be3-85c8-32735c01e66b",
    "",
    "!@#$%^&*()_+={}[]:;'<>?,./`~",
    "now valid 255 character value-at vero eosd accusamus et iusto odio dignissimos ducimus qui blanditiis praesentium voluptatum deleniti atque corrupti quos dolores et quas molestias excepturi sint occaecati cupiditate non provident, similique sunt in culpa",
  )

  val badLabelKeys = List(
    "",
    "key with characters more than 255-at vero eos et accusamus et iusto odio dignissimos ducimus qui blanditiis praesentium voluptatum deleniti atque corrupti quos dolores et quas molestias excepturi sint occaecati cupiditate non provident, similique sunt in culpa"
  )

  val badLabelValues = List(
    "value with characters more than 255-at vero eos et accusamus et iusto odio dignissimos ducimus qui blanditiis praesentium voluptatum deleniti atque corrupti quos dolores et quas molestias excepturi sint occaecati cupiditate non provident, similique sunt in culpa"
  )

  goodLabelKeys foreach { key =>
    it should s"validate a good label key '$key'" in {
      Label.validateLabelKey(key) should be(Valid(key))
    }
  }

  goodLabelValues foreach { value =>
    it should s"validate a good label value '$value'" in {
      Label.validateLabelValue(value) should be(Valid(value))
    }
  }

  badLabelKeys foreach { key =>
    it should s"not validate a bad label key $key" in {
      Label.validateLabelKey(key) match {
        case Invalid(_) => // Good!
        case Valid(_) => fail(s"Label key validation succeeded but should have failed.")
      }
    }
  }

  badLabelValues foreach { value =>
    it should s"not validate a bad label value $value" in {
      Label.validateLabelValue(value) match {
        case Invalid(_) => // Good!
        case Valid(_) => fail(s"Label value validation succeeded but should have failed.")
      }
    }
  }
}
