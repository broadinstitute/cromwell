package centaur.cwl

import io.circe.Json
import spray.json.{JsArray, JsNumber, JsString, JsValue}
import cwl.{MyriadOutputType, File => CwlFile}
import shapeless.{Inl, Poly1}
import cwl.{CwlType, MyriadOutputType, File => CwlFile}
import io.circe.Json
import io.circe.generic.auto._
import io.circe.literal._
import io.circe.refined._
import io.circe.shapes._
import io.circe.syntax._
import shapeless.{Inl, Poly1}
import spray.json.{JsNumber, JsString, JsValue}
import _root_.cwl._
import cromwell.core.path.PathBuilder

//Take cromwell's outputs and format them as expected by the spec
object OutputManipulator extends Poly1 {

  //In an Ideal world I'd return a Coproduct of these types and leave the asJson-ing to the handleOutput
  def resolveOutput(jsValue: JsValue, pathBuilder: PathBuilder, mot: MyriadOutputType): Json  = {
    mot.fold(this).apply(jsValue, pathBuilder)
  }

  private def resolveOutputViaInnerType(mot: MyriadOutputInnerType)(jsValue: JsValue, pathBuilder: PathBuilder): Json = {
    (jsValue, mot) match {
      //CWL expects quite a few enhancements to the File structure, hence...
      case (JsString(metadata), Inl(CwlType.File)) =>

        val path = pathBuilder.build(metadata).get

        CwlFile(
          location = Option(path.name),
          checksum = Option("sha1$" + path.sha1.toLowerCase),
          size = Option(path.size)
        ).asJson
      case (JsNumber(metadata), Inl(CwlType.Long)) => metadata.longValue.asJson
      case (JsNumber(metadata), Inl(CwlType.Float)) => metadata.floatValue.asJson
      case (JsNumber(metadata), Inl(CwlType.Double)) => metadata.doubleValue.asJson
      case (JsNumber(metadata), Inl(CwlType.Int)) => metadata.intValue.asJson
      case (JsString(metadata), Inl(CwlType.String)) => metadata.asJson
      case (JsArray(metadata), tpe) if tpe.select[OutputArraySchema].isDefined =>
        (for {
          schema <- tpe.select[OutputArraySchema]
          items = schema.items
          innerType <- items.select[MyriadOutputInnerType]
          outputJson = metadata.map(m => resolveOutputViaInnerType(innerType)(m, pathBuilder)).asJson
        } yield outputJson).getOrElse(throw new RuntimeException(s"We currently do not support output arrays with ${tpe.select[OutputArraySchema].get.items} inner type"))
      case (json, tpe) => throw new RuntimeException(s" we currently do not support outputs of $json and type $tpe")
    }
  }

  implicit def moit: Case.Aux[MyriadOutputInnerType, (JsValue, PathBuilder) => Json] = at[MyriadOutputInnerType] {
    resolveOutputViaInnerType(_)
  }

  implicit def amoit: Case.Aux[Array[MyriadOutputInnerType], (JsValue, PathBuilder) => Json] =
    at[Array[MyriadOutputInnerType]] {
      amoit =>
        resolveOutputViaInnerType(amoit.head)
    }
}
