package wdl4s.cwl

import io.circe.Decoder
import io.circe._


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
