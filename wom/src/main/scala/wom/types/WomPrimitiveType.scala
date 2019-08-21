package wom.types

trait WomPrimitiveType extends WomType {

  lazy val coercionMap: Map[WomType, Seq[WomType]] = Map(
    // From type -> To type
    WomAnyType -> Seq(WomStringType, WomIntegerType, WomFloatType, WomSingleFileType, WomBooleanType),
    WomStringType -> Seq(WomStringType, WomIntegerType, WomFloatType, WomSingleFileType, WomBooleanType),
    WomUnlistedDirectoryType -> Seq(WomStringType, WomUnlistedDirectoryType),
    WomSingleFileType -> Seq(WomStringType, WomSingleFileType),
    WomGlobFileType -> Seq(WomStringType, WomGlobFileType),
    WomIntegerType -> Seq(WomStringType, WomIntegerType, WomFloatType),
    WomBigDecimalType -> Seq(WomStringType, WomIntegerType, WomFloatType, WomBigDecimalType),
    WomFloatType -> Seq(WomStringType, WomFloatType),
    WomBooleanType -> Seq(WomStringType, WomBooleanType)
  )

  override def typeSpecificIsCoerceableFrom(otherType: WomType): Boolean = {
    coercionMap.get(otherType) match {
      case Some(types) => types contains this
      case None => false
    }
  }
}
