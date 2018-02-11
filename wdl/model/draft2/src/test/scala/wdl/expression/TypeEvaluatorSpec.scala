package wdl.expression

import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.{FlatSpec, Matchers}
import wdl.SampleWdl
import wdl.draft2.model.expression.WdlStandardLibraryFunctionsType
import wdl.draft2.model.types.WdlCallOutputsObjectType
import wdl.draft2.model.{WdlExpression, WdlNamespaceWithWorkflow}
import wom.types._
import wom.values.WomPrimitive

import scala.util.{Failure, Success, Try}

class TypeEvaluatorSpec extends FlatSpec with Matchers {
  val expr: String => WdlExpression = WdlExpression.fromString
  val namespace = WdlNamespaceWithWorkflow.load(SampleWdl.ThreeStep.workflowSource(), Seq.empty).get

  def noLookup(String: String): WomType = fail("No identifiers should be looked up in this test")

  def identifierLookup(String: String): WomType = {
    String match {
      case "cgrep" => WdlCallOutputsObjectType(namespace.workflow.calls.find(_.unqualifiedName == "cgrep").get)
      case "ps" => WdlCallOutputsObjectType(namespace.workflow.calls.find(_.unqualifiedName == "ps").get)
    }
  }

  def identifierEval(exprStr: String): WomPrimitiveType = expr(exprStr).evaluateType(identifierLookup, new WdlStandardLibraryFunctionsType).asInstanceOf[Try[WomPrimitiveType]].get
  def identifierEvalError(exprStr: String): Unit = {
    expr(exprStr).evaluateType(identifierLookup, new WdlStandardLibraryFunctionsType).asInstanceOf[Try[WomPrimitive]] match {
      case Failure(_) => // Expected
      case Success(badValue) => fail(s"Operation was supposed to fail, instead I got value: $badValue")
    }
  }

  private def operate(lhs: WomType, op: String, rhs: WomType): Try[WomType] = op match {
    case "+" => lhs.add(rhs)
    case "-" => lhs.subtract(rhs)
    case "*" => lhs.multiply(rhs)
    case "/" => lhs.divide(rhs)
    case "%" => lhs.mod(rhs)
    case "==" => lhs.equals(rhs)
    case "!=" => lhs.notEquals(rhs)
    case "<" => lhs.lessThan(rhs)
    case "<=" => lhs.lessThanOrEqual(rhs)
    case ">" => lhs.greaterThan(rhs)
    case ">=" => lhs.greaterThanOrEqual(rhs)
    case "||" => lhs.or(rhs)
    case "&&" => lhs.and(rhs)
    case _ => fail(s"unexpected operator: $op")
  }

  val validOperations = Table(
    ("lhs", "op", "rhs", "result"),
    (WomIntegerType, "+", WomIntegerType, WomIntegerType),
    (WomIntegerType, "+", WomFloatType, WomFloatType),
    (WomIntegerType, "+", WomStringType, WomStringType),
    (WomIntegerType, "-", WomIntegerType, WomIntegerType),
    (WomIntegerType, "-", WomFloatType, WomFloatType),
    (WomIntegerType, "*", WomIntegerType, WomIntegerType),
    (WomIntegerType, "*", WomFloatType, WomFloatType),
    (WomIntegerType, "/", WomIntegerType, WomIntegerType),
    (WomIntegerType, "/", WomFloatType, WomFloatType),
    (WomIntegerType, "%", WomIntegerType, WomIntegerType),
    (WomIntegerType, "%", WomFloatType, WomFloatType),
    (WomIntegerType, "==", WomIntegerType, WomBooleanType),
    (WomIntegerType, "==", WomFloatType, WomBooleanType),
    (WomIntegerType, "!=", WomIntegerType, WomBooleanType),
    (WomIntegerType, "!=", WomFloatType, WomBooleanType),
    (WomIntegerType, "<", WomIntegerType, WomBooleanType),
    (WomIntegerType, "<", WomFloatType, WomBooleanType),
    (WomIntegerType, "<=", WomIntegerType, WomBooleanType),
    (WomIntegerType, "<=", WomFloatType, WomBooleanType),
    (WomIntegerType, ">", WomIntegerType, WomBooleanType),
    (WomIntegerType, ">", WomFloatType, WomBooleanType),
    (WomIntegerType, ">=", WomIntegerType, WomBooleanType),
    (WomIntegerType, ">=", WomFloatType, WomBooleanType),
    (WomFloatType, "+", WomIntegerType, WomFloatType),
    (WomFloatType, "+", WomFloatType, WomFloatType),
    (WomFloatType, "+", WomStringType, WomStringType),
    (WomFloatType, "-", WomIntegerType, WomFloatType),
    (WomFloatType, "-", WomFloatType, WomFloatType),
    (WomFloatType, "*", WomIntegerType, WomFloatType),
    (WomFloatType, "*", WomFloatType, WomFloatType),
    (WomFloatType, "/", WomIntegerType, WomFloatType),
    (WomFloatType, "/", WomFloatType, WomFloatType),
    (WomFloatType, "%", WomIntegerType, WomFloatType),
    (WomFloatType, "%", WomFloatType, WomFloatType),
    (WomFloatType, "==", WomIntegerType, WomBooleanType),
    (WomFloatType, "==", WomFloatType, WomBooleanType),
    (WomFloatType, "!=", WomIntegerType, WomBooleanType),
    (WomFloatType, "!=", WomFloatType, WomBooleanType),
    (WomFloatType, "<", WomIntegerType, WomBooleanType),
    (WomFloatType, "<", WomFloatType, WomBooleanType),
    (WomFloatType, "<=", WomIntegerType, WomBooleanType),
    (WomFloatType, "<=", WomFloatType, WomBooleanType),
    (WomFloatType, ">", WomIntegerType, WomBooleanType),
    (WomFloatType, ">", WomFloatType, WomBooleanType),
    (WomFloatType, ">=", WomIntegerType, WomBooleanType),
    (WomFloatType, ">=", WomFloatType, WomBooleanType),
    (WomStringType, "+", WomIntegerType, WomStringType),
    (WomStringType, "+", WomFloatType, WomStringType),
    (WomStringType, "+", WomStringType, WomStringType),
    (WomStringType, "+", WomSingleFileType, WomStringType),
    (WomStringType, "==", WomStringType, WomBooleanType),
    (WomStringType, "!=", WomStringType, WomBooleanType),
    (WomStringType, "<", WomStringType, WomBooleanType),
    (WomStringType, "<=", WomStringType, WomBooleanType),
    (WomStringType, ">", WomStringType, WomBooleanType),
    (WomStringType, ">=", WomStringType, WomBooleanType),
    (WomSingleFileType, "+", WomStringType, WomSingleFileType),
    (WomSingleFileType, "==", WomStringType, WomBooleanType),
    (WomSingleFileType, "==", WomSingleFileType, WomBooleanType),
    (WomSingleFileType, "!=", WomStringType, WomBooleanType),
    (WomSingleFileType, "!=", WomSingleFileType, WomBooleanType),
    (WomSingleFileType, "<=", WomStringType, WomBooleanType),
    (WomSingleFileType, "<=", WomSingleFileType, WomBooleanType),
    (WomSingleFileType, ">=", WomStringType, WomBooleanType),
    (WomSingleFileType, ">=", WomSingleFileType, WomBooleanType),
    (WomBooleanType, "==", WomBooleanType, WomBooleanType),
    (WomBooleanType, "!=", WomBooleanType, WomBooleanType),
    (WomBooleanType, "<", WomBooleanType, WomBooleanType),
    (WomBooleanType, "<=", WomBooleanType, WomBooleanType),
    (WomBooleanType, ">", WomBooleanType, WomBooleanType),
    (WomBooleanType, ">=", WomBooleanType, WomBooleanType),
    (WomBooleanType, "||", WomBooleanType, WomBooleanType),
    (WomBooleanType, "&&", WomBooleanType, WomBooleanType)
  )

  val invalidOperations = Table(
    ("lhs", "op", "rhs"),
    (WomIntegerType, "+", WomSingleFileType),
    (WomIntegerType, "+", WomBooleanType),
    (WomIntegerType, "-", WomStringType),
    (WomIntegerType, "-", WomSingleFileType),
    (WomIntegerType, "-", WomBooleanType),
    (WomIntegerType, "*", WomStringType),
    (WomIntegerType, "*", WomSingleFileType),
    (WomIntegerType, "*", WomBooleanType),
    (WomIntegerType, "/", WomStringType),
    (WomIntegerType, "/", WomSingleFileType),
    (WomIntegerType, "/", WomBooleanType),
    (WomIntegerType, "%", WomStringType),
    (WomIntegerType, "%", WomSingleFileType),
    (WomIntegerType, "%", WomBooleanType),
    (WomIntegerType, "==", WomStringType),
    (WomIntegerType, "==", WomSingleFileType),
    (WomIntegerType, "==", WomBooleanType),
    (WomIntegerType, "!=", WomStringType),
    (WomIntegerType, "!=", WomSingleFileType),
    (WomIntegerType, "!=", WomBooleanType),
    (WomIntegerType, "<", WomStringType),
    (WomIntegerType, "<", WomSingleFileType),
    (WomIntegerType, "<", WomBooleanType),
    (WomIntegerType, "<=", WomStringType),
    (WomIntegerType, "<=", WomSingleFileType),
    (WomIntegerType, "<=", WomBooleanType),
    (WomIntegerType, ">", WomStringType),
    (WomIntegerType, ">", WomSingleFileType),
    (WomIntegerType, ">", WomBooleanType),
    (WomIntegerType, ">=", WomStringType),
    (WomIntegerType, ">=", WomSingleFileType),
    (WomIntegerType, ">=", WomBooleanType),
    (WomIntegerType, "||", WomIntegerType),
    (WomIntegerType, "||", WomFloatType),
    (WomIntegerType, "||", WomStringType),
    (WomIntegerType, "||", WomSingleFileType),
    (WomIntegerType, "||", WomBooleanType),
    (WomIntegerType, "&&", WomIntegerType),
    (WomIntegerType, "&&", WomFloatType),
    (WomIntegerType, "&&", WomStringType),
    (WomIntegerType, "&&", WomSingleFileType),
    (WomIntegerType, "&&", WomBooleanType),
    (WomFloatType, "+", WomSingleFileType),
    (WomFloatType, "+", WomBooleanType),
    (WomFloatType, "-", WomStringType),
    (WomFloatType, "-", WomSingleFileType),
    (WomFloatType, "-", WomBooleanType),
    (WomFloatType, "*", WomStringType),
    (WomFloatType, "*", WomSingleFileType),
    (WomFloatType, "*", WomBooleanType),
    (WomFloatType, "/", WomStringType),
    (WomFloatType, "/", WomSingleFileType),
    (WomFloatType, "/", WomBooleanType),
    (WomFloatType, "%", WomStringType),
    (WomFloatType, "%", WomSingleFileType),
    (WomFloatType, "%", WomBooleanType),
    (WomFloatType, "==", WomStringType),
    (WomFloatType, "==", WomSingleFileType),
    (WomFloatType, "==", WomBooleanType),
    (WomFloatType, "!=", WomStringType),
    (WomFloatType, "!=", WomSingleFileType),
    (WomFloatType, "!=", WomBooleanType),
    (WomFloatType, "<", WomStringType),
    (WomFloatType, "<", WomSingleFileType),
    (WomFloatType, "<", WomBooleanType),
    (WomFloatType, "<=", WomStringType),
    (WomFloatType, "<=", WomSingleFileType),
    (WomFloatType, "<=", WomBooleanType),
    (WomFloatType, ">", WomStringType),
    (WomFloatType, ">", WomSingleFileType),
    (WomFloatType, ">", WomBooleanType),
    (WomFloatType, ">=", WomStringType),
    (WomFloatType, ">=", WomSingleFileType),
    (WomFloatType, ">=", WomBooleanType),
    (WomFloatType, "||", WomIntegerType),
    (WomFloatType, "||", WomFloatType),
    (WomFloatType, "||", WomStringType),
    (WomFloatType, "||", WomSingleFileType),
    (WomFloatType, "||", WomBooleanType),
    (WomFloatType, "&&", WomIntegerType),
    (WomFloatType, "&&", WomFloatType),
    (WomFloatType, "&&", WomStringType),
    (WomFloatType, "&&", WomSingleFileType),
    (WomFloatType, "&&", WomBooleanType),
    (WomStringType, "+", WomBooleanType),
    (WomStringType, "-", WomIntegerType),
    (WomStringType, "-", WomFloatType),
    (WomStringType, "-", WomStringType),
    (WomStringType, "-", WomSingleFileType),
    (WomStringType, "-", WomBooleanType),
    (WomStringType, "*", WomIntegerType),
    (WomStringType, "*", WomFloatType),
    (WomStringType, "*", WomStringType),
    (WomStringType, "*", WomSingleFileType),
    (WomStringType, "*", WomBooleanType),
    (WomStringType, "/", WomIntegerType),
    (WomStringType, "/", WomFloatType),
    (WomStringType, "/", WomStringType),
    (WomStringType, "/", WomSingleFileType),
    (WomStringType, "/", WomBooleanType),
    (WomStringType, "%", WomIntegerType),
    (WomStringType, "%", WomFloatType),
    (WomStringType, "%", WomStringType),
    (WomStringType, "%", WomSingleFileType),
    (WomStringType, "%", WomBooleanType),
    (WomStringType, "==", WomIntegerType),
    (WomStringType, "==", WomFloatType),
    (WomStringType, "==", WomSingleFileType),
    (WomStringType, "==", WomBooleanType),
    (WomStringType, "!=", WomIntegerType),
    (WomStringType, "!=", WomFloatType),
    (WomStringType, "!=", WomSingleFileType),
    (WomStringType, "!=", WomBooleanType),
    (WomStringType, "<", WomIntegerType),
    (WomStringType, "<", WomFloatType),
    (WomStringType, "<", WomSingleFileType),
    (WomStringType, "<", WomBooleanType),
    (WomStringType, "<=", WomIntegerType),
    (WomStringType, "<=", WomFloatType),
    (WomStringType, "<=", WomSingleFileType),
    (WomStringType, "<=", WomBooleanType),
    (WomStringType, ">", WomIntegerType),
    (WomStringType, ">", WomFloatType),
    (WomStringType, ">", WomSingleFileType),
    (WomStringType, ">", WomBooleanType),
    (WomStringType, ">=", WomIntegerType),
    (WomStringType, ">=", WomFloatType),
    (WomStringType, ">=", WomSingleFileType),
    (WomStringType, ">=", WomBooleanType),
    (WomStringType, "||", WomIntegerType),
    (WomStringType, "||", WomFloatType),
    (WomStringType, "||", WomStringType),
    (WomStringType, "||", WomSingleFileType),
    (WomStringType, "||", WomBooleanType),
    (WomStringType, "&&", WomIntegerType),
    (WomStringType, "&&", WomFloatType),
    (WomStringType, "&&", WomStringType),
    (WomStringType, "&&", WomSingleFileType),
    (WomStringType, "&&", WomBooleanType),
    (WomSingleFileType, "+", WomSingleFileType),
    (WomSingleFileType, "+", WomIntegerType),
    (WomSingleFileType, "+", WomFloatType),
    (WomSingleFileType, "+", WomBooleanType),
    (WomSingleFileType, "-", WomIntegerType),
    (WomSingleFileType, "-", WomFloatType),
    (WomSingleFileType, "-", WomStringType),
    (WomSingleFileType, "-", WomSingleFileType),
    (WomSingleFileType, "-", WomBooleanType),
    (WomSingleFileType, "*", WomIntegerType),
    (WomSingleFileType, "*", WomFloatType),
    (WomSingleFileType, "*", WomStringType),
    (WomSingleFileType, "*", WomSingleFileType),
    (WomSingleFileType, "*", WomBooleanType),
    (WomSingleFileType, "/", WomIntegerType),
    (WomSingleFileType, "/", WomFloatType),
    (WomSingleFileType, "/", WomStringType),
    (WomSingleFileType, "/", WomSingleFileType),
    (WomSingleFileType, "/", WomBooleanType),
    (WomSingleFileType, "%", WomIntegerType),
    (WomSingleFileType, "%", WomFloatType),
    (WomSingleFileType, "%", WomStringType),
    (WomSingleFileType, "%", WomSingleFileType),
    (WomSingleFileType, "%", WomBooleanType),
    (WomSingleFileType, "==", WomIntegerType),
    (WomSingleFileType, "==", WomFloatType),
    (WomSingleFileType, "==", WomBooleanType),
    (WomSingleFileType, "!=", WomIntegerType),
    (WomSingleFileType, "!=", WomFloatType),
    (WomSingleFileType, "!=", WomBooleanType),
    (WomSingleFileType, "<", WomIntegerType),
    (WomSingleFileType, "<", WomFloatType),
    (WomSingleFileType, "<", WomStringType),
    (WomSingleFileType, "<", WomSingleFileType),
    (WomSingleFileType, "<", WomBooleanType),
    (WomSingleFileType, "<=", WomIntegerType),
    (WomSingleFileType, "<=", WomFloatType),
    (WomSingleFileType, "<=", WomBooleanType),
    (WomSingleFileType, ">", WomIntegerType),
    (WomSingleFileType, ">", WomFloatType),
    (WomSingleFileType, ">", WomStringType),
    (WomSingleFileType, ">", WomSingleFileType),
    (WomSingleFileType, ">", WomBooleanType),
    (WomSingleFileType, ">=", WomIntegerType),
    (WomSingleFileType, ">=", WomFloatType),
    (WomSingleFileType, ">=", WomBooleanType),
    (WomSingleFileType, "||", WomIntegerType),
    (WomSingleFileType, "||", WomFloatType),
    (WomSingleFileType, "||", WomStringType),
    (WomSingleFileType, "||", WomSingleFileType),
    (WomSingleFileType, "||", WomBooleanType),
    (WomSingleFileType, "&&", WomIntegerType),
    (WomSingleFileType, "&&", WomFloatType),
    (WomSingleFileType, "&&", WomStringType),
    (WomSingleFileType, "&&", WomSingleFileType),
    (WomSingleFileType, "&&", WomBooleanType),
    (WomBooleanType, "+", WomIntegerType),
    (WomBooleanType, "+", WomFloatType),
    (WomBooleanType, "+", WomStringType),
    (WomBooleanType, "+", WomSingleFileType),
    (WomBooleanType, "+", WomBooleanType),
    (WomBooleanType, "-", WomIntegerType),
    (WomBooleanType, "-", WomFloatType),
    (WomBooleanType, "-", WomStringType),
    (WomBooleanType, "-", WomSingleFileType),
    (WomBooleanType, "-", WomBooleanType),
    (WomBooleanType, "*", WomIntegerType),
    (WomBooleanType, "*", WomFloatType),
    (WomBooleanType, "*", WomStringType),
    (WomBooleanType, "*", WomSingleFileType),
    (WomBooleanType, "*", WomBooleanType),
    (WomBooleanType, "/", WomIntegerType),
    (WomBooleanType, "/", WomFloatType),
    (WomBooleanType, "/", WomStringType),
    (WomBooleanType, "/", WomSingleFileType),
    (WomBooleanType, "/", WomBooleanType),
    (WomBooleanType, "%", WomIntegerType),
    (WomBooleanType, "%", WomFloatType),
    (WomBooleanType, "%", WomStringType),
    (WomBooleanType, "%", WomSingleFileType),
    (WomBooleanType, "%", WomBooleanType),
    (WomBooleanType, "==", WomIntegerType),
    (WomBooleanType, "==", WomFloatType),
    (WomBooleanType, "==", WomStringType),
    (WomBooleanType, "==", WomSingleFileType),
    (WomBooleanType, "!=", WomIntegerType),
    (WomBooleanType, "!=", WomFloatType),
    (WomBooleanType, "!=", WomStringType),
    (WomBooleanType, "!=", WomSingleFileType),
    (WomBooleanType, "<", WomIntegerType),
    (WomBooleanType, "<", WomFloatType),
    (WomBooleanType, "<", WomStringType),
    (WomBooleanType, "<", WomSingleFileType),
    (WomBooleanType, "<=", WomIntegerType),
    (WomBooleanType, "<=", WomFloatType),
    (WomBooleanType, "<=", WomStringType),
    (WomBooleanType, "<=", WomSingleFileType),
    (WomBooleanType, ">", WomIntegerType),
    (WomBooleanType, ">", WomFloatType),
    (WomBooleanType, ">", WomStringType),
    (WomBooleanType, ">", WomSingleFileType),
    (WomBooleanType, ">=", WomIntegerType),
    (WomBooleanType, ">=", WomFloatType),
    (WomBooleanType, ">=", WomStringType),
    (WomBooleanType, ">=", WomSingleFileType),
    (WomBooleanType, "||", WomIntegerType),
    (WomBooleanType, "||", WomFloatType),
    (WomBooleanType, "||", WomStringType),
    (WomBooleanType, "||", WomSingleFileType),
    (WomBooleanType, "&&", WomIntegerType),
    (WomBooleanType, "&&", WomFloatType),
    (WomBooleanType, "&&", WomStringType),
    (WomBooleanType, "&&", WomSingleFileType)
  )

  forAll (validOperations) { (lhs, op, rhs, expectedType) =>
    it should s"validate the output type for the expression: $lhs $op $rhs = $expectedType" in {
      operate(lhs, op, rhs) shouldEqual Success(expectedType)
    }
  }

  forAll (invalidOperations) { (lhs, op, rhs) =>
    it should s"not allow the expression: $lhs $op $rhs" in {
      operate(lhs, op, rhs) should be(a[Failure[_]])
    }
  }

  "Expression Evaluator with Object as LHS" should "Lookup object string attribute" in {
    identifierEval("cgrep.count") shouldEqual WomIntegerType
  }
  it should "Lookup object integer attribute" in {
    identifierEval("ps.procs") shouldEqual WomSingleFileType
  }
  it should "Error if key doesn't exist" in {
    identifierEvalError("ps.badkey")
  }
}
