package wdl4s

import wdl4s.AstTools.{AstNodeName, EnhancedAstNode}
import wdl4s.WdlExpression.ScopedLookupFunction
import wdl4s.expression.{NoFunctions, WdlFunctions, WdlStandardLibraryFunctionsType}
import wdl4s.parser.MemoryUnit
import wdl4s.parser.WdlParser.{Ast, AstList, SyntaxError}
import wdl4s.types.{WdlIntegerType, WdlStringType}
import wdl4s.util.TryUtil
import wdl4s.values._

import scala.collection.JavaConverters._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

case class RuntimeAttributes(attrs: Map[String, WdlExpression]) {
  import RuntimeAttributes._

  def evaluate(lookup: ScopedLookupFunction, functions: WdlFunctions[WdlValue]): Map[String, Try[WdlValue]] = {
    attrs map { case (k, v) =>
      k -> v.evaluate(lookup, functions).flatMap(validateRuntimeValue(k, _))
    }
  }

  private def validateRuntimeValue(key: String, value: WdlValue): Try[WdlValue] = {
    key match {
      case "memory" => validateMemoryValue(value).map(m => WdlString(m.toString))
      case _ => Success(value)
    }
  }
}

object RuntimeAttributes {
  def apply(ast: Ast, declarations: Seq[Declaration]): RuntimeAttributes = {
    val asts = ast.findAsts(AstNodeName.Runtime)
    if (asts.size > 1) throw new UnsupportedOperationException("Only one runtime block may be defined per task")
    val kvPairAsts = asts.headOption.map(_.getAttribute("map").asInstanceOf[AstList].asScala.toVector.map(_.asInstanceOf[Ast]))
    val runtimeAttributeMap = kvPairAsts match {
      case Some(vector) => vector.map(ast => processRuntimeAttribute(ast, declarations)).toMap
      case None => Map.empty[String, Try[WdlExpression]]
    }

    // .get below to throw the exception if validation didn't pass
    TryUtil.sequenceMap(runtimeAttributeMap).map(RuntimeAttributes(_)).get
  }

  def validateMemoryValue(value: WdlValue): Try[MemorySize] = {
    value match {
      case i: WdlInteger => Success(MemorySize(i.value.toDouble, MemoryUnit.Bytes))
      case s: WdlString => MemorySize.parse(s.valueString).map(m => m.to(MemoryUnit.GB))
      case other => Failure(new UnsupportedOperationException("Valid memory values are either strings (e.g. '8 GB') or integers"))
    }
  }

  private def processRuntimeAttribute(ast: Ast, declarations: Seq[Declaration]): (String, Try[WdlExpression]) = {
    val key = ast.getAttribute("key").sourceString
    val expression = WdlExpression(ast.getAttribute("value"))
    val validatedExpression = key match {
      case "memory" => validateMemoryExpression(expression, declarations)
      case _ => Success(expression)
    }
    key -> validatedExpression
  }

  private def validateMemoryExpression(expression: WdlExpression, declarations: Seq[Declaration]): Try[WdlExpression] = {
    // .get below because lookup functions should throw exceptions if they could not lookup the variable
    def lookupType(n: String) = declarations.find(_.name == n).map(_.wdlType).get
    expression.evaluateType(lookupType, new WdlStandardLibraryFunctionsType) match {
      case Success(wdlType) if !Seq(WdlIntegerType, WdlStringType).contains(wdlType) =>
        Failure(new SyntaxError("Expecting the 'memory' runtime attribute to be either a String or an Int"))
      case Success(wdlType) => validateStaticMemoryExpression(expression)
      case Failure(ex) => Failure(ex)
    }
  }

  /** If an expression can be evaluated without any lookup function or WDL function implementations,
    * Then validate that the static value is properly formatted.  If the expression cannot be statically
    * evaluated at compile time, then defer the validation until runtime (i.e. return Success).
    *
    * For example, the expression `"8" + "GB"` can be statically evaluated to `"8GB"`, which can be
    * determined to be correctly formatted at compile time.  However, the expression `x + "GB` cannot
    * be statically evaluated because we may not have a value for `x` yet, so return Success(expression)
    * and check the value when the expression can be evaluated.
    */
  private def validateStaticMemoryExpression(expression: WdlExpression): Try[WdlExpression] = {
    expression.evaluate(NoLookup, NoFunctions) match {
      case Success(value) => validateMemoryValue(value).map(_ => expression)
      case Failure(ex) => Success(expression)
    }
  }
}
