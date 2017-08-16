package cromwell.webservice

import org.scalatest.{FlatSpec, Matchers}
import spray.json._
import spray.json.DefaultJsonProtocol._

class PartialWorkflowSourcesSpec extends FlatSpec with Matchers {
  it should "succesfully merge and override multiple input files" in {
    val input1 = Map("wf.a1" -> "hello", "wf.a2" -> "world").toJson.toString
    val input2 = Map.empty[String, String].toJson.toString
    val overrideInput1 = Map("wf.a2" -> "universe").toJson.toString
    val allInputs = PartialWorkflowSources.mergeMaps(Seq(Option(input1), Option(input2), Option(overrideInput1)))

    allInputs.fields.keys should contain allOf("wf.a1", "wf.a2")
    allInputs.fields("wf.a2") should be(JsString("universe"))
  }
}
