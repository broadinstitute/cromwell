package wdl.command

import wdl.AstTools.EnhancedAstNode
import wdl._
import wdl.exception.VariableNotFoundException
import wdl.expression.WdlFunctions
import wdl4s.parser.WdlParser.{Ast, SyntaxError, Terminal}
import wom.OptionalNotSuppliedException
import wom.types.WdlOptionalType
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

  private def withPostProcesser(f: WdlValue => WdlValue, valueFunctions: WdlFunctions[WdlValue]) = new WdlFunctions[WdlValue] {
    override def getFunction(name: String): WdlFunction = {
      val originalFunction = valueFunctions.getFunction(name)
      (seq: Seq[Try[WdlValue]]) => originalFunction(seq) map f
    }
  }
}

case class ParameterCommandPart(attributes: Map[String, String], expression: WdlExpression) extends WdlCommandPart {
  def attributesToString: String = if (attributes.nonEmpty) attributes.map({case (k,v) => s"$k=${WdlString(v).toWdlString}"}).mkString(" ") + " " else ""
  override def toString: String = "${" + s"$attributesToString${expression.toWdlString}" + "}"

  override def instantiate(declarations: Seq[Declaration], inputs: Map[String, WdlValue], functions: WdlFunctions[WdlValue], valueMapper: (WdlValue) => WdlValue): String = {
    // This is a safety net.
    // In Cromwell's production code, optional declarations are always passed to instantiate, as WdlOptionalValue.none(type) if necessary.
    def lookupDeclaration(s: String) = declarations.collectFirst {
      // The backtick syntax (`s`) allows us to equality-check 's' against the match/case result:
      case Declaration(wdlType: WdlOptionalType, `s`, _, _, _) => wdlType.none
    } getOrElse { throw VariableNotFoundException(s) }

    val lookup: String => WdlValue = (s: String) => valueMapper(inputs.getOrElse(s, lookupDeclaration(s)))

    val value = expression.evaluate(lookup, ParameterCommandPart.withPostProcesser(valueMapper, functions)) match {
      case Success(v) => v match {
        case WdlOptionalValue(_, opt) => opt.getOrElse(defaultString)
        case _ => v
      }
      case Failure(OptionalNotSuppliedException(_)) => defaultString
      case Failure(f) => throw new UnsupportedOperationException(s"Could not evaluate expression: ${expression.toWdlString}", f)
    }

    valueMapper(value) match {
      case b: WdlBoolean if attributes.contains("true") && attributes.contains("false") => if (b.value) attributes.get("true").head else attributes.get("false").head
      case p: WdlPrimitive => p.valueString
      case a: WdlArray if attributes.contains("sep") => a.value.map(_.valueString).mkString(attributes.get("sep").head)
      case _: WdlArray => throw new UnsupportedOperationException(s"Expression '${expression.toString}' evaluated to an Array but no 'sep' was specified")
      case _ => throw new UnsupportedOperationException(s"Could not string-ify value: $value")
    }
  }

  private def defaultString = {
    if (attributes.contains("default")) {
      WdlString(attributes.get("default").head)
    } else {
      WdlString("")
    }
  }
}
