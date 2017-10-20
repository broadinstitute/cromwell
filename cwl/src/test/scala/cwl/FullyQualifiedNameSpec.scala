package cwl

import org.scalatest.{FlatSpec, Matchers}

class FullyQualifiedNameSpec extends FlatSpec with Matchers {

  "workflow input id " should "get filename and id" in {
    val wfid = WorkflowInputId("file#id")

    wfid.fileName shouldBe "file"
    wfid.inputId shouldBe "id"
  }

  "workflow step output id " should "get filename, id, and step" in {
    val wfid = WorkflowStepInputOrOutputId("file#step/id")

    wfid.fileName shouldBe "file"
    wfid.ioId shouldBe "id"
    wfid.stepId shouldBe "step"
  }
}
