package wdl4s.cwl

import io.circe.Decoder
import io.circe._
import io.circe.generic.auto._

import io.circe.shapes._
import io.circe.generic.auto._
import shapeless.Coproduct
import cats.syntax.either._
import eu.timepit.refined.string._
import eu.timepit.refined._
import io.circe.refined._
import io.circe.literal._
import io.circe.syntax._
import cats.data.ValidatedNel
import cats.data.Validated._
import cats.syntax.traverse._
import cats.instances.list._
import cats.syntax.option._
import lenthall.validation.ErrorOr.ErrorOr

object Implicits {

  implicit val cwlTypeDecoder = Decoder.enumDecoder(CwlType)
  implicit val cwlVersionDecoder = Decoder.enumDecoder(CwlVersion)
  implicit val scatterMethodDecoder = Decoder.enumDecoder(ScatterMethod)
  implicit val linkMergeMethodDecoder = Decoder.enumDecoder(LinkMergeMethod)

  implicit val cwlTypeEncoder = Encoder.enumEncoder(CwlType)
  implicit val cwlVersionEncoder = Encoder.enumEncoder(CwlVersion)
  implicit val scatterMethodEncoder = Encoder.enumEncoder(ScatterMethod)
  implicit val linkMergeMethodEncoder = Encoder.enumEncoder(LinkMergeMethod)

  implicit def enumerationEncoder[V <: Enumeration#Value]: Encoder[V] = (value: V) => Json.fromString(value.toString)
}
