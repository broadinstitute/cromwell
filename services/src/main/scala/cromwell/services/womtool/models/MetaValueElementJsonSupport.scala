package cromwell.services.womtool.models

import cats.syntax.functor._
//import io.circe.{Encoder, Json}
import io.circe.{ Decoder, Encoder }, io.circe.generic.auto._
import io.circe.syntax._

import wom.callable.MetaValueElement
import wom.callable.MetaValueElement._

object MetaValueElementJsonSupport {

  implicit val metaValueElementEncoder: Encoder[MetaValueElement] = Encoder.instance {
    case MetaValueElementNull => throw new Exception("null not supported yet")
    case a @ MetaValueElementBoolean(_) => a.asJson
    case a @ MetaValueElementFloat(_) => a.asJson
    case a @ MetaValueElementInteger(_) => a.asJson
    case a @ MetaValueElementString(_) => a.asJson
    case MetaValueElementObject(_) => throw new Exception("object not supported yet")
    case MetaValueElementArray(_) => throw new Exception("array not supported yet")
  }

  implicit val metaValueElementDecoder: Decoder[MetaValueElement] =
    List[Decoder[MetaValueElement]](
//      Decoder[MetaValueElementNull].widen,
      Decoder[MetaValueElementBoolean].widen,
      Decoder[MetaValueElementFloat].widen,
      Decoder[MetaValueElementInteger].widen,
      Decoder[MetaValueElementString].widen,
    ).reduceLeft(_ or _)
}
