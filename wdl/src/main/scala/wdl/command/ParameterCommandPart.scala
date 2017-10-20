package wdl.command

import wdl.AstTools.EnhancedAstNode
import wdl._
import wdl.exception.VariableNotFoundException
import wdl.expression.WdlFunctions
import wdl4s.parser.WdlParser.{Ast, SyntaxError, Terminal}
import wom.OptionalNotSuppliedException
import wom.types.WomOptionalType
import wom.values._

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object ParameterCommandPart {
  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): ParameterCommandPart = {
    val attributes = ast.getAttribute("attributes").astListAsVector map { a =>
      val ast = a.asInstanceOf[Ast]
      (ast.getAttribute("key").sourceString, ast.getAttribute("value").sourceString)
    } toMap
    val expression = WdlExpression(ast.getAttribute("expr"))
    if ((attributes.contains("true") && !attributes.contains("false")) || (attributes.contains("false") && !attributes.contains("true"))) {
      // .head because we can't get here without there being at least one attribute
      val firstAttr = ast.getAttribute("attributes").astListAsVector.head.asInstanceOf[Ast].getAttribute("key").asInstanceOf[Terminal]
      throw new SyntaxError(wdlSyntaxErrorFormatter.trueAndFalseAttributesAreRequired(firstAttr))
    }
    new ParameterCommandPart(attributes, expression)
  }

  private def withPostProcesser(f: WomValue => WomValue, valueFunctions: WdlFunctions[WomValue]) = new WdlFunctions[WomValue] {
    override def getFunction(name: String): WdlFunction = {
      val originalFunction = valueFunctions.getFunction(name)
      (seq: Seq[Try[WomValue]]) => originalFunction(seq) map f
    }
  }
}

case class ParameterCommandPart(attributes: Map[String, String], expression: WdlExpression) extends WdlCommandPart {
  def attributesToString: String = if (attributes.nonEmpty) attributes.map({case (k,v) => s"$k=${WomString(v).toWomString}"}).mkString(" ") + " " else ""
  override def toString: String = "${" + s"$attributesToString${expression.toWomString}" + "}"

  override def instantiate(declarations: Seq[Declaration], inputs: Map[String, WomValue], functions: WdlFunctions[WomValue], valueMapper: (WomValue) => WomValue): String = {
    // This is a safety net.
    // In Cromwell's production code, optional declarations are always passed to instantiate, as WdlOptionalValue.none(type) if necessary.
    def lookupDeclaration(s: String) = declarations.collectFirst {
      // The backtick syntax (`s`) allows us to equality-check 's' against the match/case result:
      case Declaration(womType: WomOptionalType, `s`, _, _, _) => womType.none
    } getOrElse { throw VariableNotFoundException(s) }

    val lookup: String => WomValue = (s: String) => valueMapper(inputs.getOrElse(s, lookupDeclaration(s)))

    val value = expression.evaluate(lookup, ParameterCommandPart.withPostProcesser(valueMapper, functions)) match {
      case Success(v) => v match {
        case WomOptionalValue(_, opt) => opt.getOrElse(defaultString)
        case _ => v
      }
      case Failure(OptionalNotSuppliedException(_)) => defaultString
      case Failure(f) => throw new UnsupportedOperationException(s"Could not evaluate expression: ${expression.toWomString}", f)
    }

    valueMapper(value) match {
      case b: WomBoolean if attributes.contains("true") && attributes.contains("false") => if (b.value) attributes.get("true").head else attributes.get("false").head
      case p: WomPrimitive => p.valueString
      case a: WomArray if attributes.contains("sep") => a.value.map(_.valueString).mkString(attributes.get("sep").head)
      case _: WomArray => throw new UnsupportedOperationException(s"Expression '${expression.toString}' evaluated to an Array but no 'sep' was specified")
      case _ => throw new UnsupportedOperationException(s"Could not string-ify value: $value")
    }
  }

  private def defaultString = {
    if (attributes.contains("default")) {
      WomString(attributes.get("default").head)
    } else {
      WomString("")
    }
  }
}
