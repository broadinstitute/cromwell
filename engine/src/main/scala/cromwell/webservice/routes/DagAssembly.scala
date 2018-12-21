package cromwell.webservice.routes

import cats.implicits._
import cromwell.core.WorkflowId
import spray.json._

object DagAssembly extends DefaultJsonProtocol {

  final case class RgbAverage(rTotal: Long, rCount: Int, gTotal: Long, gCount: Int, bTotal: Long, bCount: Int) {
    private def hexify(total: Long, count: Int): String = {
      val result = (total / count).toHexString.toUpperCase
      if (result.length == 0) "00" + result else if (result.length == 1) "0" + result else result
    }
    def dotString = s"#${hexify(rTotal, rCount)}${hexify(gTotal, gCount)}${hexify(bTotal, bCount)}"
  }
  object RgbAverage {
    def apply(r: Int, g: Int, b: Int): RgbAverage = RgbAverage(r.longValue(), 1, g.longValue(), 1, b.longValue(), 1)
  }
  implicit val rgbAverageMonoid = cats.derived.MkMonoid[RgbAverage]

  def dagFromMetadata(workflowId: WorkflowId, metadata: JsObject): String = {
    val noDagFound = s"""|digraph noDAG {
                         |  #rankdir=LR;
                         |  compound=true;
                         |
                         |  subgraph cluster_0 {
                         |  graph [color=red label="Cannot find DAG in metadata" penwidth=3]
                         |  node [shape=rectangle color=red]
                         |
                         |  A [label="NO graph information found in the metadata of this workflow ID ${workflowId.toString}" fontname="times bold"]
                         |  B [label="But the workflow ID *is* valid"]
                         |
                         |  C [label="Perhaps this workflow hasn't been started yet? Cromwell needs to analyze"]
                         |  D [label="Cromwell has a rate limit on how fast it picks up new workflows"]
                         |  E [label="Parsing the workflow can take a few seconds - especially if there are a lot of remote imports"]
                         |
                         |  F [label="Perhaps the workflow was invalid?"]
                         |  G [label="If Cromwell cannot construct an interal representation, it cannot produce the graph view"]
                         |
                         |  H [label="Perhaps this workflow was started before Cromwell added the 'dag' endpoint?"]
                         |  I [label="Cromwell now writes a 'dag' field to metadata as soon as it parses a workflow"]
                         |  J [label="If Cromwell never wrote the data, we cannot display the graph diagram"]
                         |  }
                         |
                         |  A -> B
                         |  B -> C
                         |  B -> F
                         |  B -> H
                         |  C -> D
                         |  C -> E
                         |  F -> G
                         |  H -> I
                         |  H -> J
                         |}""".stripMargin

    if (metadata.fields.contains("dag")) metadata.fields("dag").convertTo[String] else noDagFound
  }


  /**
    * Convert the 'calls' object from metadata into a mapping of call name to status color.
    */
  def callColorsFromMetadata(metadata: JsObject): Map[String, JsString] = {

    val calls = metadata.fields("calls").asJsObject.fields

    calls map { case (fqn, value) =>
      val name = fqn.split('.').last

      val array = value.convertTo[List[JsValue]]
      val fillColor = if (array.nonEmpty) {
        array.foldMap { callEntry =>
          callEntry.asJsObject.fields("executionStatus").convertTo[String] match {
            case "Done" => RgbAverage(0x00, 0x80, 0x01)
            case "QueuedInCromwell" => RgbAverage(0x7F, 0xFF, 0xD4)
            case "Running" => RgbAverage(0x64, 0x95, 0xED)
            case _ => RgbAverage(0xFF, 0x00, 0x00)
          }
        }
      } else {
        RgbAverage(0xD3, 0xD3, 0xD3)
      }

      s"CALL_$name" -> JsString(fillColor.dotString)
    }
  }
}
