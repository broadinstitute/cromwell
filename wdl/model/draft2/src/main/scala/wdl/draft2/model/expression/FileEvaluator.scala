package wdl.draft2.model.expression

import common.util.TryUtil
import wdl.draft2.model.AstTools.EnhancedAstNode
import wdl.draft2.model.TernaryIf
import wdl.draft2.model.WdlExpression._
import wdl.draft2.parser.WdlParser.{Ast, AstNode}
import wdl.shared.model.expression.FileEvaluatorUtil
import wom.WomExpressionException
import wom.types._
import wom.values._

import scala.util.{Failure, Success, Try}

class FileEvaluatorWdlFunctions

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
case class FileEvaluator(valueEvaluator: ValueEvaluator, coerceTo: WomType = WomAnyType) {

  private def evalValue(ast: AstNode): Try[WomValue] = valueEvaluator.evaluate(ast)

  private def evalValueToWdlFile(ast: AstNode): Try[WomFile] = {
    evalValue(ast) match {
      case Success(p: WomPrimitive) => Success(WomSingleFile(p.valueString))
      case Success(_) => Failure(new WomExpressionException(s"Expecting a primitive type from AST:\n${ast.toPrettyString}"))
      case Failure(e) => Failure(e)
    }
  }

  def evaluate(ast: AstNode, anticipatedType: WomType = coerceTo): Try[Seq[WomFile]] = {
    valueEvaluator.evaluate(ast) match {
        // If the ast can successfully be evaluated, then just find the WdlFiles that it contains
      case Success(v) => Success(FileEvaluatorUtil.findFilesToDelocalize(v, coerceTo))
      case Failure(_) => evaluateRecursive(ast, anticipatedType)
    }
  }

  /**
    * Recursively traverse the ast and collect asts that evaluate successfully to a WdlFile
    * or that take a WdlFile as a parameter, and that parameter can be evaluated and is a WdlPrimitive
    * For example, assuming that the command produces a file called "out", take the following output section
    *   output {
    *     String x = read_string("out")
    *   }
    *
    * Purely from the output signature (String x) there is no way to know that this task is expected to produce a file.
    * If we try to evaluate read_string("out") before the task is run it will most likely fail as "out" hasn't been created yet.
    * So we look at what type of ast read_string is, in this case it's a "FunctionCallWithOneFileParameter",
    * that means that whatever the parameter is it should evaluates to a file.
    * We then try to evaluate the "out" string which directly produces WdlString("out").
    * Because we know that read_string takes a file parameter, WdlString("out") is transformed to a WdlFile (see evalValueToWdlFile)
    * We can now deduce that this task is expected to produce a WdlFile called "out"
    */
  private def evaluateRecursive(ast: AstNode, anticipatedType: WomType): Try[Seq[WomFile]] = {
    ast match {
      case a: Ast if a.isGlobFunctionCall =>
        evalValueToWdlFile(a.params.head) map { wdlFile => Seq(WomGlobFile(wdlFile.value)) }
      case a: Ast if a.isFunctionCallWithFirstParameterBeingFile =>
        evalValueToWdlFile(a.params.head) match {
          case Success(v) => Success(Seq(v))
          case _ => evaluateRecursive(a.params.head, anticipatedType)
        }
      case a: Ast if a.isFunctionCall =>
        TryUtil.sequence(a.params map { p => evaluateRecursive(p, anticipatedType) }) map { _.flatten }
      case a: Ast if a.isBinaryOperator =>
        evalValueToWdlFile(a) match {
          case Success(f: WomFile) => Success(Seq(f))
          case _ => for {
            left <- evaluateRecursive(a.getAttribute("lhs"), anticipatedType)
            right <- evaluateRecursive(a.getAttribute("rhs"), anticipatedType)
          } yield left ++ right
        }
      case a: Ast if a.isUnaryOperator =>
        evalValueToWdlFile(a) match {
          case Success(f: WomFile) => Success(Seq(f))
          case _ =>
            evaluateRecursive(a.getAttribute("expression"), anticipatedType) match {
              case Success(a:Seq[WomFile]) => Success(a)
              case _ => Failure(new WomExpressionException(s"Could not evaluate:\n${a.toPrettyString}"))
            }
        }
      case TernaryIf(condition, ifTrue, ifFalse) =>
        for {
          filesInCondition <- evaluate(condition, WomBooleanType)
          filesIfTrue <- evaluate(ifTrue)
          filesIfFalse <- evaluate(ifFalse)
        } yield filesInCondition ++ filesIfTrue ++ filesIfFalse
      case a: Ast if a.isArrayOrMapLookup =>
        evalValue(a) match {
          case Success(f: WomFile) => Success(Seq(f))
          case _ => for {
            left <- evaluateRecursive(a.getAttribute("lhs"), anticipatedType)
            right <- evaluateRecursive(a.getAttribute("rhs"), anticipatedType)
          } yield left ++ right
        }
      case a: Ast if a.isMemberAccess =>
        evalValue(a) match {
          case Success(f: WomFile) => Success(Seq(f))
          case _ => Success(Seq.empty[WomFile])
        }
      case a: Ast if a.isArrayLiteral =>
        anticipatedType match {
          case WomArrayType(memberType) =>
            val values = a.getAttribute("values").astListAsVector.map(x => evaluate(x, memberType))
            TryUtil.sequence(values) match {
              case Success(v) => Success(v.flatten)
              case f => f.map(_.flatten)
            }
          case _ => Failure(new WomExpressionException(s"Failed to parse $a for files. Found an unexpected Array literal but anticipated a ${anticipatedType.stableName}"))
        }
      case a: Ast if a.isTupleLiteral =>
        val unevaluatedElements = a.getAttribute("values").astListAsVector
        if (unevaluatedElements.size == 1) {
          evaluateRecursive(unevaluatedElements.head, anticipatedType)
        } else if (unevaluatedElements.size == 2) {
          anticipatedType match {
            case WomPairType(leftType, rightType) =>
              for {
                left <- evaluate(unevaluatedElements.head, leftType)
                right <- evaluate(unevaluatedElements(1), rightType)
              } yield left ++ right
            case _ => Failure(new WomExpressionException(s"Failed to parse $a for files. Found an unexpected Pair literal but anticipated a ${anticipatedType.stableName}"))
          }

        } else {
          Failure(new WomExpressionException(s"WDL does not currently support tuples with n > 2: ${a.toPrettyString}"))
        }
      case a: Ast if a.isMapLiteral =>
        anticipatedType match {
          case WomMapType(keyType, valueType) =>
            val evaluatedMap = a.getAttribute("map").astListAsVector map { kv =>
              val key = evaluate(kv.asInstanceOf[Ast].getAttribute("key"), keyType)
              val value = evaluate(kv.asInstanceOf[Ast].getAttribute("value"), valueType)
              key -> value
            }

            val flattenedTries = evaluatedMap flatMap { case (k,v) => Seq(k,v) }
            flattenedTries partition {_.isSuccess} match {
              case (_, failures) if failures.nonEmpty =>
                val message = failures.collect { case f: Failure[_] => f.exception.getMessage }.mkString("\n")
                Failure(new WomExpressionException(s"Could not evaluate expression:\n$message"))
              case (successes, _) =>
                Success(successes.flatMap(_.get))
            }
          case _ => Failure(new WomExpressionException(s"Failed to parse $a for files. Found an unexpected Map literal but anticipated a ${anticipatedType.stableName}"))
        }


      case _ => Success(Seq.empty[WomFile])
    }
  }
}

