package wdl4s.wdl.types

abstract class WdlPrimitiveType extends WdlType {
  lazy val coercionMap: Map[WdlType, Seq[WdlType]] = Map(
    // From type -> To type
    WdlStringType -> Seq(WdlStringType, WdlIntegerType, WdlFloatType, WdlFileType, WdlBooleanType),
    WdlFileType -> Seq(WdlStringType, WdlFileType),
    WdlIntegerType -> Seq(WdlStringType, WdlIntegerType, WdlFloatType),
    WdlFloatType -> Seq(WdlStringType, WdlFloatType),
    WdlBooleanType -> Seq(WdlStringType, WdlBooleanType)
  )

  override def isCoerceableFrom(otherType: WdlType): Boolean = {
    coercionMap.get(otherType) match {
      case Some(types) => types contains this
      case None => false
    }
  }
}

