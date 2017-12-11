package centaur.cwl

import io.circe.Json
import spray.json.{JsNumber, JsString, JsValue}
import _root_.cwl._
import cwl.{MyriadOutputType, File => CwlFile}
import cats.effect.IO
import shapeless.{Inl, Poly1}
import better.files.File
import cwl.{CwlType, MyriadOutputType, File => CwlFile}
import io.circe.Json
import io.circe.generic.auto._
import io.circe.literal._
import io.circe.refined._
import io.circe.shapes._
import io.circe.syntax._
import shapeless.{Inl, Poly1}
import spray.json.{JsNumber, JsString, JsValue}
import fs2.io.file.readAll
import fs2.hash.sha1
import _root_.cwl._
import better.files.File

//Take cromwell's outputs and format them as expected by the spec
object OutputManipulator extends Poly1{

  //In an Ideal world I'd return a Coproduct of these types and leave the asJson-ing to the handleOutput
  def resolveOutput(jsValue: JsValue, mot: MyriadOutputType): Json  = {
    mot.fold(this).apply(jsValue)
  }

  private def resolveOutputViaInnerType(mot: MyriadOutputInnerType)(jsValue:JsValue) : Json  = {
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

  implicit def moit: Case.Aux[MyriadOutputInnerType, JsValue => Json] = at[MyriadOutputInnerType] {
    resolveOutputViaInnerType(_)
  }

  implicit def amoit: Case.Aux[Array[MyriadOutputInnerType], JsValue => Json] =
    at[Array[MyriadOutputInnerType]] {
      amoit =>
        resolveOutputViaInnerType(amoit.head)
    }

}
