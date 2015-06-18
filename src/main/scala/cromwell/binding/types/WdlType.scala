package cromwell.binding.types

import cromwell.binding.values.WdlValue

import scala.util.{Failure, Try}

trait WdlType {

  /**
   * Method to be overridden by implementation classes defining a partial function
   * for the conversion of raw input values to specific implementation class value types.
   * i.e.  `WdlBooleanType` should define a partial function that knows how to
   * construct `WdlBoolean`s for inputs of supported types and contents.  Values for which
   * the partial function is not defined are assumed to not be convertible to the target type.
   */
  protected def coercion: PartialFunction[Any, WdlValue]

  /**
   * Public interface for a `Try`-wrapped conversion of an input of type `Any` to
   * a `WdlValue`.
   */
  def coerceRawValue(any: Any): Try[WdlValue] = {
    if (!coercion.isDefinedAt(any)) {
      Failure(new IllegalArgumentException(s"No coercion defined from '$any' of type '${any.getClass}' to ${getClass.getSimpleName}."))
    } else {
      Try(coercion(any))
    }
  }

  def toWdlString: String
}

object WdlType {

  private lazy val wdlTypes = Seq(
    WdlBooleanType, WdlExpressionType, WdlFileType, WdlFloatType,
    WdlIntegerType, WdlNamespaceType, WdlObjectType, WdlStringType
  )

  def fromWdlString(wdlString: String): WdlType = wdlTypes.find(_.toWdlString == wdlString).getOrElse(
    throw new NoSuchElementException(s"No such WdlType: $wdlString"))
}
