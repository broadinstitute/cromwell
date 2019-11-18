package cromwell.api.model

import org.scalatest.{FlatSpec, Matchers}

class WaasDescriptionJsonSupportSpec extends FlatSpec with Matchers {

  it should "deserialize invalid result JSON" in {
    val resultJson =
      """{
        |  "valid": false,
        |  "errors": [
        |    "Failed to import workflow sub_workflow_aborted_import.wdl.:\\nBad import sub_workflow_aborted_import.wdl: Failed to resolve 'sub_workflow_aborted_import.wdl' using resolver: 'http importer (no 'relative-to' origin)' (reason 1 of 1): Relative path"
        |  ],
        |  "validWorkflow": false,
        |  "name": "",
        |  "inputs": [],
        |  "outputs": [],
        |  "images": [],
        |  "submittedDescriptorType": {},
        |  "importedDescriptorTypes": [],
        |  "meta": {},
        |  "parameterMeta": {},
        |  "isRunnableWorkflow": false
        |}""".stripMargin


    import cromwell.api.model.WorkflowDescriptionJsonSupport._
    import spray.json._

    val jsonAst = resultJson.parseJson
    val deserialized = jsonAst.convertTo[WaasDescription]

    deserialized.valid should be(false)
    deserialized.errors should be(List("""Failed to import workflow sub_workflow_aborted_import.wdl.:\nBad import sub_workflow_aborted_import.wdl: Failed to resolve 'sub_workflow_aborted_import.wdl' using resolver: 'http importer (no 'relative-to' origin)' (reason 1 of 1): Relative path"""))
    deserialized.validWorkflow should be(false)
    deserialized.name should be("")
    deserialized.inputs should be(List.empty)
    deserialized.outputs should be(List.empty)
    deserialized.images should be(List.empty)
    deserialized.submittedDescriptorType should be(WaasWorkflowDescriptorType(None, None))
    deserialized.importedDescriptorTypes should be(List.empty)
    deserialized.meta should be(JsObject.empty)
    deserialized.parameterMeta should be(JsObject.empty)
    deserialized.isRunnableWorkflow should be(false)
  }
}
