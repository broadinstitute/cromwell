package centaur.cwl

import cats.effect.IO
import centaur.api.CentaurCromwellClient
import cromwell.api.model.SubmittedWorkflow
import cwl.{CwlDecoder, CwlType, MyriadOutputType, File => CwlFile}
import io.circe.Json
import io.circe.generic.auto._
import io.circe.literal._
import io.circe.refined._
import io.circe.shapes._
import io.circe.syntax._
import shapeless.{Inl, Inr, Poly1}
import spray.json.{JsNumber, JsObject, JsString, JsValue}
import fs2.io.file.readAll
import fs2.hash.sha1
import _root_.cwl._
import better.files.File

import scalaz.syntax.std.map._

object Outputs {

  //When the string returned is not valid JSON, it is effectively an exception as CWL runner expects JSON to be returned
  def handleOutput(submittedWorkflow: SubmittedWorkflow): String = {
    val metadata: Map[String, JsValue] = CentaurCromwellClient.metadata(submittedWorkflow).get.value

    //Sorry for all the nesting, but spray json JsValue doesn't have optional access methods like Argonaut/circe,
    //thus no for comprehensions for us :(
    metadata.get("submittedFiles.workflow") match {
      case Some(JsString(workflow)) =>

        val cwl = CwlDecoder.decodeTopLevelCwl(workflow)

        cwl.value.attempt.unsafeRunSync() match {
          case Right(Right(cwl)) =>

            CentaurCromwellClient.outputs(submittedWorkflow).get.outputs match {
              case JsObject(map) =>
                val typeMap: Map[String, MyriadOutputType] = cwl.fold(CwlOutputsFold)
                val mungeTypeMap = typeMap.mapKeys(stripTypeMapKey)

                val mungeOutputMap = map.mapKeys(stripOutputKey)

                mungeOutputMap.
                  //This lets us operate on the values of the output values and types for a particular output key
                  intersectWith(mungeTypeMap)(resolveOutput).
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

  //Ids come out of SALAD pre-processing with a filename prepended.  This gets rid of it
  def stripTypeMapKey(key: String): String = key.substring(key.lastIndexOf("#") + 1, key.length)

  //Ids come out of Cromwell with a prefix, separated by a ".".  This takes everything to the right,
  //as CWL wants it
  def stripOutputKey(key: String): String = key.substring(key.lastIndexOf(".") + 1, key.length)

  //In an Ideal world I'd return a Coproduct of these types and leave the asJson-ing to the handleOutput
  def resolveOutput(jsValue: JsValue, mot: MyriadOutputType): Json  = {
    mot match {
      case Inl(moit) => resolveOutputViaInnerType(jsValue, moit)
      case Inr(Inl(arrayMoit)) => resolveOutputViaInnerType(jsValue, arrayMoit.head)
      case _ => ??? //Impossible!
    }
  }
  def resolveOutputViaInnerType(jsValue: JsValue, mot: MyriadOutputInnerType): Json  = {
    (jsValue, mot) match {
      //CWL expects quite a few enhancements to the File structure, hence...
      case (JsString(metadata), Inl(CwlType.File)) =>

        val fileName = File(metadata)

        val file = readAll[IO](fileName.path, 65536)
        val hash = sha1[IO]

        val bytes = (file through hash).runLog.unsafeRunSync

        //hexify the bytes into a string
        //credit: https://stackoverflow.com/questions/2756166/what-is-are-the-scala-ways-to-implement-this-java-byte-to-hex-class
        val hex = bytes.map {
          b => String.format("%02X", new Integer(b & 0xff))
        }.mkString.toLowerCase

        CwlFile(
          location = Option(fileName.name),
          checksum = Option("sha1$" + hex),
          size = Option(fileName.size)
        ).asJson
      case (JsNumber(metadata), Inl(CwlType.Long)) => metadata.longValue.asJson
      case (JsNumber(metadata), Inl(CwlType.Float)) => metadata.floatValue.asJson
      case (JsNumber(metadata), Inl(CwlType.Double)) => metadata.doubleValue.asJson
      case (JsNumber(metadata), Inl(CwlType.Int)) => metadata.intValue.asJson
      case (JsString(metadata), Inl(CwlType.String)) => metadata.asJson
      case (json, tpe) => throw new RuntimeException(s" we currently do not support outputs of $json and type $tpe")
    }
  }
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
    _.outputs.map(output => output.id -> output.`type`).toMap
  }
}
