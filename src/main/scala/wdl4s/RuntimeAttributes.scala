package wdl4s

import wdl4s.AstTools.{AstNodeName, EnhancedAstNode}
import wdl4s.expression.NoFunctions
import wdl4s.values._
import wdl4s.parser.WdlParser.{Ast, AstList}
import wdl4s.parser.MemorySize

import scala.collection.JavaConverters._
import scala.language.postfixOps
import scalaz.Scalaz._
import scalaz._

case class RuntimeAttributes(attrs: Map[String, Seq[String]])

object RuntimeAttributes {
  def apply(ast: Ast): RuntimeAttributes = {
    val asts = ast.findAsts(AstNodeName.Runtime)
    if (asts.size > 1) throw new UnsupportedOperationException("Only one runtime block may be defined per task")
    val astList = asts.headOption map { _.getAttribute("map").asInstanceOf[AstList] }
    val attrMap = astList map processRuntimeAttributes getOrElse Map.empty[String, Seq[String]]
    attrMap.get("memory") flatMap { _.headOption } foreach validateMemory 

    RuntimeAttributes(attrMap)
  }

  /**
    * Checks that a provided memory value is valid, if not throws an IllegalArgumentException
    */
  private def validateMemory(value: String): Unit = {
    validateMemStringInGb (value) match {
      case Failure(f) => throw new IllegalArgumentException("Invalid memory value specified:\n" + f.toList.mkString("\n"))
      case Success(_) => ()
    }
  }

  private def validateMemStringInGb(mem: String): ValidationNel[String, Double] = {
    val memoryPattern = """(\d+)\s*(\w+)""".r
    mem match {
      case memoryPattern(amountString, unitString) =>
        val amount = amountString.parseLong leftMap { _.getMessage } toValidationNel
        val unit = validateMemoryUnit(unitString)
        (amount |@| unit) { (a, u) => MemorySize.GB.fromBytes(u.toBytes(a)) }
      case _ => s"$mem should be of the form X Unit where X is a number, e.g. 8 GB".failureNel
    }
  }

  private def validateMemoryUnit(unit: String): ValidationNel[String, MemorySize] = {
    MemorySize.values find { _.suffixes.contains(unit) } match {
      case Some(s) => s.successNel
      case None => s"$unit is an invalid memory unit".failureNel
    }
  }

  private def processRuntimeAttributes(astList: AstList): Map[String, Seq[String]] = {
    astList.asScala.toVector map { a => processRuntimeAttribute(a.asInstanceOf[Ast]) } toMap
  }

  private def processRuntimeAttribute(ast: Ast): (String, Seq[String]) = {
    val key = ast.getAttribute("key").sourceString
    val seq = Option(ast.getAttribute("value")) map { valAttr =>
      WdlExpression.evaluate(valAttr, NoLookup, NoFunctions) match {
        case scala.util.Success(wdlPrimitive: WdlPrimitive) => Seq(wdlPrimitive.valueString)
        case scala.util.Success(wdlArray: WdlArray) => wdlArray.value.map(_.valueString)
        case scala.util.Success(wdlValue: WdlValue) =>
          throw new IllegalArgumentException(
            s"WdlType not supported for $key, ${wdlValue.wdlType}: ${wdlValue.valueString}")
        case null =>
          throw new IllegalArgumentException(s"value was null: $key}")
        case scala.util.Failure(f) => throw f
      }
    } getOrElse Seq.empty

    (key, seq)
  }
}
