package cromwell.binding.expression

import cromwell.binding.WdlExpression._
import cromwell.binding.WdlExpressionException
import cromwell.binding.types._
import cromwell.binding.values._
import cromwell.parser.WdlParser.{Ast, AstNode}
import cromwell.binding.AstTools.EnhancedAstNode
import cromwell.util.TryUtil

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

class FileEvaluatorWdlFunctions extends WdlFunctions[Seq[WdlFile]]

/**
 * This evaluator will take a WdlExpression and determine all of the static files that
 * are referenced in this expression.
 *
 * This utilizes a ValueEvaluator to try to statically evaluate parts of the expression
 * into WdlValues.  See the evalValueToWdlFile() function.
 *
 * The coerceTo parameter is for the case where an output might be coerceable to a File.
 * Consider the following output section for a task:
 *
 * output {
 *   File a = "a.txt"
 *   String b = "b.txt"
 * }
 *
 * In the first case, the coerceTo would be WdlFileType and the calling evaluate("a.txt") would
 * return a Seq(WdlFile("a.txt")).  In the second case, coerceTo would be WdlStringType and
 * evaluate("b.txt") would return Seq().
 */
case class FileEvaluator(valueEvaluator: ValueEvaluator, coerceTo: WdlType = WdlAnyType) extends Evaluator {
  override type T = Seq[WdlFile]
  override def lookup = (s:String) => Seq.empty[WdlFile]
  override val functions = new FileEvaluatorWdlFunctions()

  private def evalValue(ast: AstNode): Try[WdlValue] = valueEvaluator.evaluate(ast)

  private def evalValueToWdlFile(ast: AstNode): Try[Seq[WdlFile]] = {
    evalValue(ast) match {
      case Success(p: WdlPrimitive) => Success(Seq(WdlFile(p.valueString)))
      case Success(_) => Failure(new WdlExpressionException(s"Expecting a primitive type from AST:\n${ast.toPrettyString}"))
      case Failure(e) => Failure(e)
    }
  }

  private def findWdlFiles(value: WdlValue, coerce: Boolean = true): Seq[WdlFile] = {
    val coercedValue = if (coerce) coerceTo.coerceRawValue(value) else Success(value)
    coercedValue match {
      case Success(f: WdlFile) => Seq(f)
      case Success(a: WdlArray) =>
        a.value.flatMap(findWdlFiles(_, coerce=false))
      case Success(m: WdlMap) =>
        m.value flatMap { case (k, v) => Seq(k, v) } flatMap(findWdlFiles(_, coerce=false)) toSeq
      case _ => Seq.empty[WdlFile]
    }
  }

  override def evaluate(ast: AstNode): Try[Seq[WdlFile]] = {
    /**
     * First check if the top-level expression evaluates to a
     * literal value.  If it does, return all WdlFiles referenced
     * in that literal value
     */
    valueEvaluator.evaluate(ast) match {
      case Success(v) => Success(findWdlFiles(v))
      case Failure(ex) => evaluateRecursive(ast)
    }
  }

  /**
   * The pattern below is to try to evaluate all ASTs to a static value first, via
   * evalValueToWdlFile().  If this succeeds, return that value.  Otherwise, call
   * this function recursively.
   */
  private def evaluateRecursive(ast: AstNode): Try[Seq[WdlFile]] = {
    ast match {
      case a: Ast if a.isFunctionCallWithOneFileParameter =>
        evalValueToWdlFile(a.params.head) match {
          case Success(v) => Success(v)
          case _ => evaluateRecursive(a.params.head)
        }
      case a: Ast if a.isBinaryOperator =>
        evalValueToWdlFile(a) match {
          case Success(f: WdlFile) => Success(Seq(f))
          case _ => (evaluateRecursive(a.getAttribute("lhs")), evaluateRecursive(a.getAttribute("rhs"))) match {
            case (Success(a:Seq[WdlFile]), Success(b:Seq[WdlFile])) => Success(a ++ b)
            case _ => Failure(new WdlExpressionException(s"Could not evaluate:\n${a.toPrettyString}"))
          }
        }
      case a: Ast if a.isUnaryOperator =>
        evalValueToWdlFile(a) match {
          case Success(f: WdlFile) => Success(Seq(f))
          case _ =>
            evaluateRecursive(a.getAttribute("expression")) match {
              case Success(a:Seq[WdlFile]) => Success(a)
              case _ => Failure(new WdlExpressionException(s"Could not evaluate:\n${a.toPrettyString}"))
            }
        }
      case a: Ast if a.isArrayOrMapLookup =>
        evalValue(a) match {
          case Success(f: WdlFile) => Success(Seq(f))
          case _ => evaluateRecursive(a.getAttribute("rhs"))
        }
      case a: Ast if a.isMemberAccess =>
        evalValue(a) match {
          case Success(f: WdlFile) => Success(Seq(f))
          case _ => Success(Seq.empty[WdlFile])
        }
      case a: Ast if a.isArrayLiteral =>
        val values = a.getAttribute("values").astListAsVector.map(evaluateRecursive)
        TryUtil.sequence(values) match {
          case Success(v) => Success(v.flatten)
          case f => f.map(_.flatten)
        }
      case a: Ast if a.isMapLiteral =>
        val evaluatedMap = a.getAttribute("map").astListAsVector map { kv =>
          val key = evaluateRecursive(kv.asInstanceOf[Ast].getAttribute("key"))
          val value = evaluateRecursive(kv.asInstanceOf[Ast].getAttribute("value"))
          key -> value
        }

        val flattenedTries = evaluatedMap flatMap { case (k,v) => Seq(k,v) }
        flattenedTries partition {_.isSuccess} match {
          case (_, failures) if failures.nonEmpty =>
            val message = failures.collect { case f: Failure[_] => f.exception.getMessage }.mkString("\n")
            Failure(new WdlExpressionException(s"Could not evaluate expression:\n$message"))
          case (successes, _) =>
            Success(successes.flatMap(_.get))
        }
      case _ => Success(Seq.empty[WdlFile])
    }
  }
}

