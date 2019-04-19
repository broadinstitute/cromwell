package cromwell.backend.google.pipelines.common

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}

class GoogleLabelsSpec extends FlatSpec with Matchers {

  behavior of "GoogleLabels"

  /**
    * In the format 'to validate', 'expected result'
    */
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

  googleLabelConversions foreach { case (label: String, conversion: String) =>
    it should s"not validate the bad label key '$label'" in {
      GoogleLabels.validateLabelRegex(label) match {
        case Invalid(_) => // Good!
        case Valid(_) => fail(s"Label validation succeeded but should have failed.")
      }
    }

    it should s"convert the bad label string '$label' into the safe label string '$conversion'" in {
      GoogleLabels.safeGoogleName(label) should be(conversion)
    }
  }
}
