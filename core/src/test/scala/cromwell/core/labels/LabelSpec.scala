package cromwell.core.labels

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}

class LabelSpec extends FlatSpec with Matchers {

  behavior of "Labels"

  /**
    * In the format 'to validate', 'expected result'
    */
  val goodLabelKeys = List(
    "cromwell-root-workflow-id",
    "cromwell-11f2468c-39d6-4be3-85c8-32735c01e66b",
    "just-the-right-length-just-the-right-length-just-the-right-leng"
  )

  val goodLabelValues = List(
    "11f2468c-39d6-4be3-85c8-32735c01e66b",
    ""
  )

  val badLabelKeys = List(
    "11f2468c-39d6-4be3-85c8-32735c01e66b",
    "0-cromwell-root-workflow-id",
    "",
    "cromwell-root-workflow-id-",
    "0-cromwell-root-workflow-id-",
    "Cromwell-root-workflow-id"
  )

  val googleLabelConversions = List(
    "11f2468c-39d6-4be3-85c8-32735c01e66b" -> "x--11f2468c-39d6-4be3-85c8-32735c01e66b",
    "0-cromwell-root-workflow-id" -> "x--0-cromwell-root-workflow-id",
    "cromwell-root-workflow-id-" -> "cromwell-root-workflow-id---x",
    "0-cromwell-root-workflow-id-" -> "x--0-cromwell-root-workflow-id---x",
    "Cromwell-root-workflow-id" -> "cromwell-root-workflow-id",
    "cromwell_root_workflow_id" -> "cromwell-root-workflow-id",
    "too-long-too-long-too-long-too-long-too-long-too-long-too-long-t" -> "too-long-too-long-too-long-too---g-too-long-too-long-too-long-t",
    "0-too-long-and-invalid-too-long-and-invalid-too-long-and-invali+" -> "x--0-too-long-and-invalid-too----nvalid-too-long-and-invali---x"
  )

  goodLabelKeys foreach { label =>
    it should s"validate a good label key '$label'" in {
      Label.validateLabelKey(label) should be(Valid(label))
    }
  }

  goodLabelValues foreach { label =>
    it should s"validate a good label value '$label'" in {
      Label.validateLabelValue(label) should be(Valid(label))
    }
  }

  googleLabelConversions foreach { case (label: String, conversion: String) =>
    it should s"not validate the bad label string '$label'" in {
      Label.validateLabelRegex(label, Label.GoogleLabelRegexPattern.r) match {
        case Invalid(_) => // Good!
        case Valid(_) => fail(s"Label validation succeeded but should have failed.")
      }
    }

    it should s"convert the bad label string '$label' into the safe label string '$conversion'" in {
      Label.safeGoogleName(label) should be(conversion)
    }
  }
}
