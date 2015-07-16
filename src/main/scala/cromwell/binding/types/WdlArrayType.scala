package cromwell.binding.types

import cromwell.binding.{WdlExpressionException, WdlFunctions, WdlExpression}
import cromwell.binding.values.{WdlInteger, WdlArray}
import cromwell.parser.WdlParser.{AstList, Ast}
import spray.json.JsArray
import scala.collection.JavaConverters._
import scala.util.{Success, Failure}

case class WdlArrayType(memberType: WdlType) extends WdlType {
  val toWdlString: String = s"Array[${memberType.toWdlString}]"

  private def coerceIterable(s: Seq[Any]): WdlArray = {
    val coerced = s.map {memberType.coerceRawValue(_).get}
    WdlArray(WdlArrayType(coerced.head.wdlType), coerced)
  }

  override protected def coercion = {
    case s: Seq[Any] if s.nonEmpty => coerceIterable(s)
    case js: JsArray if js.elements.nonEmpty => coerceIterable(js.elements)
  }
}
