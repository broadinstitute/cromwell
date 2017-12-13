package wom.types

abstract class WomPrimitiveType extends WomType {
  lazy val coercionMap: Map[WomType, Seq[WomType]] = Map(
    // From type -> To type
    WomStringType -> Seq(WomStringType, WomIntegerType, WomFloatType, WomSingleFileType, WomBooleanType),
    WomSingleFileType -> Seq(WomStringType, WomSingleFileType),
    WomIntegerType -> Seq(WomStringType, WomIntegerType, WomFloatType),
    WomFloatType -> Seq(WomStringType, WomFloatType),
    WomBooleanType -> Seq(WomStringType, WomBooleanType)
  )

  override def isCoerceableFrom(otherType: WomType): Boolean = {
    coercionMap.get(otherType) match {
      case Some(types) => types contains this
      case None => false
    }
  }
}

