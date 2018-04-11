package cwl.internal

import cwl.CommandLineTool.CommandInputParameter
import cwl.CwlType
import cwl.SchemaDefRequirement.SchemaDefTypes
import org.scalatest.{FlatSpec, Matchers}
import shapeless.Coproduct
import cwl._

class CommandPartSortingAlgorithmSpec extends FlatSpec with Matchers {

  /*
  "inputbindingsCommandPart" should "handle optional values correctly" in {
    val miit = Coproduct[cwl.MyriadInputInnerType](CwlType.Int)
    val miit2 = Coproduct[cwl.MyriadInputInnerType](CwlType.Null)


    val mit = Coproduct[MyriadInputType](Array(miit, miit2))

    val cip = CommandInputParameter(id = "x", `type` = Some(mit))

    val inputRecordField =

    val sdt = Coproduct[SchemaDefTypes](InputRecordSchema(id = "Stage", fields = Some(Array())))

    val sdr = SchemaDefRequirement(types = Array)

    CommandPartSortingAlgorithm.inputBindingsCommandPart(cip).run
  }
  */

}
