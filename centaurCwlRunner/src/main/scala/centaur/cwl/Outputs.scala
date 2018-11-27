package centaur.cwl

import centaur.api.CentaurCromwellClient
import centaur.test.metadata.WorkflowFlatMetadata._
import common.validation.IOChecked._
import cromwell.api.model.SubmittedWorkflow
import cromwell.core.path.PathBuilder
import cwl.ontology.Schema
import cwl.{Cwl, CwlDecoder, MyriadOutputType}
import io.circe.Json
import io.circe.syntax._
import scalaz.syntax.std.map._
import shapeless.Poly1
import spray.json.{JsObject, JsString, JsValue}

object Outputs {

  //When the string returned is not valid JSON, it is effectively an exception as CWL runner expects JSON to be returned
  def handleOutput(submittedWorkflow: SubmittedWorkflow, pathBuilder: PathBuilder): String = {
    val metadata: Map[String, JsValue] = CentaurCromwellClient.metadata(submittedWorkflow).unsafeRunSync().asFlat.value

    // Wrapper function to provide the right signature for `intersectWith` below.
    def outputResolver(schemaOption: Option[Schema])(jsValue: JsValue, mot: MyriadOutputType): Json = {
      OutputManipulator.resolveOutput(jsValue, pathBuilder, mot, schemaOption)
    }

    //Sorry for all the nesting, but spray json JsValue doesn't have optional access methods like Argonaut/circe,
    //thus no for comprehensions for us :(
    metadata.get("submittedFiles.workflow") match {
      case Some(JsString(workflow)) =>

        val parseCwl: IOChecked[Cwl] = CwlDecoder.decodeCwlString(
          workflow,
          submittedWorkflow.workflow.zippedImports,
          submittedWorkflow.workflow.workflowRoot
        )

        parseCwl.value.attempt.unsafeRunSync() match {
          case Right(Right(cwl)) =>

            CentaurCromwellClient.outputs(submittedWorkflow).unsafeRunSync().outputs match {
              case JsObject(map) =>
                val typeMap: Map[String, MyriadOutputType] = cwl.fold(CwlOutputsFold)
                val mungeTypeMap = typeMap.mapKeys(stripTypeMapKey)

                val mungeOutputMap = map.mapKeys(stripOutputKey)

                mungeOutputMap.
                  //This lets us operate on the values of the output values and types for a particular output key
                  intersectWith(mungeTypeMap)(outputResolver(cwl.schemaOption)).
                  //converting the whole response to Json using Circe's auto-encoder derivation
                  asJson.
                  //drop null values so that we don't print when Option == None
                  pretty(io.circe.Printer.spaces2.copy(dropNullValues = true))
              case other => s"it seems cromwell is not returning outputs as a Jsobject but is instead a $other"
            }
          case Right(Left(error)) => s"couldn't parse workflow: $workflow failed with error: $error"
          case Left(error) => s"Exception when trying to read workflow: $workflow failed with error: $error"
        }
      case Some(other) => s"received the value $other when the workflow string was expected"
      case None => "the workflow is no longer in the metadata payload, it's a problem"
    }
  }

  //Ids come out of SALAD pre-processing with a filename prepended. This gets rid of it
  // Also gets rid of the parent workflow name if present
  def stripTypeMapKey(key: String): String = key
    .substring(key.lastIndexOf("#") + 1, key.length)
    .split("/").last

  //Ids come out of Cromwell with a prefix, separated by a ".". This takes everything to the right,
  //as CWL wants it
  def stripOutputKey(key: String): String = key.substring(key.lastIndexOf(".") + 1, key.length)

}

object CwlOutputsFold extends Poly1 {
  import cwl._

  implicit def wf: Case.Aux[cwl.Workflow, Map[String, MyriadOutputType]] = at[cwl.Workflow] {
    _.outputs.map(output => output.id -> output.`type`.get).toMap
  }

  implicit def clt: Case.Aux[cwl.CommandLineTool, Map[String, MyriadOutputType]] = at[cwl.CommandLineTool] {
    _.outputs.map(output => output.id -> output.`type`.get).toMap
  }

  implicit def et: Case.Aux[cwl.ExpressionTool, Map[String, MyriadOutputType]] = at[cwl.ExpressionTool] {
    _.outputs.map(output => output.id -> output.`type`.get).toMap
  }
}

