package wdl.draft3.transforms.wdlom2wom

import common.Checked
import common.validation.Checked._
import wdl.model.draft3.elements.{PrimitiveTypeElement, TypeElement}
import wom.types.WomType

object TypeElementToType {
  def convert(a: TypeElementToWomTypeParameters): Checked[WomType] = a.typeElement match {
    case PrimitiveTypeElement(p) => p.validNelCheck
    case other => s"No implemented conversion to WOM type for '$other'".invalidNelCheck
  }

  final case class TypeElementToWomTypeParameters(typeElement: TypeElement) // availableStructs: Map[String, WomObjectType]
}


