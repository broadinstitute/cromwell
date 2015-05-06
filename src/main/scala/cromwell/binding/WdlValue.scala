package cromwell.binding

import cromwell.binding.types.{WdlFloatType, WdlIntegerType, WdlStringType, WdlType}

case class WdlValue(value: Any, wdlType: WdlType) {
  wdlType.checkCompatible(value)

  def add(rhsValue: WdlValue): WdlValue = {
    val rhs = rhsValue.value
    val lhs = value
    (wdlType, rhsValue.wdlType) match {
      case (WdlStringType, WdlIntegerType) =>
        new WdlValue(lhs.asInstanceOf[String] + rhs.asInstanceOf[Integer], WdlStringType)
      case (WdlIntegerType, WdlStringType) =>
        new WdlValue(lhs.asInstanceOf[Integer] + rhs.asInstanceOf[String], WdlStringType)
      case (WdlStringType, WdlStringType) =>
        new WdlValue(lhs.asInstanceOf[String] + rhs.asInstanceOf[String], WdlStringType)
      case (WdlIntegerType, WdlIntegerType) =>
        new WdlValue(lhs.asInstanceOf[Integer] + rhs.asInstanceOf[Integer], WdlIntegerType)
      case (WdlFloatType, WdlFloatType) =>
        new WdlValue(lhs.asInstanceOf[Float] + rhs.asInstanceOf[Float], WdlFloatType)
    }
  }
}
