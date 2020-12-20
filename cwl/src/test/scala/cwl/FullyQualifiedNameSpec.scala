package cwl

import common.assertion.CromwellTimeoutSpec
import cwl.command.ParentName
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class FullyQualifiedNameSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  implicit val parentName = ParentName.empty
  
  "workflow input id " should "get filename and id" in {
    val wfid = FileAndId("file#id")

    wfid.fileName shouldBe "file"
    wfid.id shouldBe "id"
  }

  "workflow step output id " should "get filename, id, and step" in {
    val wfid = FileStepAndId("file#step/id")

    wfid.fileName shouldBe "file"
    wfid.id shouldBe "id"
    wfid.stepId shouldBe "step"
  }
}
