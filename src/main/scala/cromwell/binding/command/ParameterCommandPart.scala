package cromwell.binding.command

import cromwell.binding.AstTools.EnhancedAstNode
import cromwell.binding._
import cromwell.binding.expression.{NoFunctions, WdlFunctions}
import cromwell.binding.values.{WdlArray, WdlPrimitive, WdlString, WdlValue}
import cromwell.parser.WdlParser.Ast

import scala.util.{Failure, Success}

object ParameterCommandPart {
  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): ParameterCommandPart = {
    val attributes = ast.getAttribute("attributes").astListAsVector map { a =>
      val ast = a.asInstanceOf[Ast]
      (ast.getAttribute("key").sourceString, ast.getAttribute("value").sourceString)
    }
    val expression = WdlExpression(ast.getAttribute("expr"))
    new ParameterCommandPart(attributes.toMap, expression)
  }
}

case class ParameterCommandPart(attributes: Map[String, String], expression: WdlExpression) extends CommandPart {
  def attributesToString: String = if (attributes.nonEmpty) attributes.map({case (k,v) => s"$k=${WdlString(v).toWdlString}"}).mkString(", ") + " " else ""
  override def toString: String = "${" + s"$attributesToString${expression.toWdlString}" + "}"

  override def instantiate(declarations: Seq[Declaration], parameters: Map[String, WdlValue], functions: WdlFunctions[WdlValue] = new NoFunctions): String = {
    val value = expression.evaluate(WdlExpression.standardLookupFunction(parameters, declarations, functions), functions) match {
      case Success(v) => v
      case Failure(f) => f match {
        case v: VariableNotFoundException => declarations.find(_.name == v.variable) match {
          /* Allow an expression to fail evaluation if one of the variables that it requires is optional (the type has ? after it, e.g. String?) */
          case Some(d) if d.postfixQuantifier.contains("?") => if (attributes.contains("default")) WdlString(attributes.get("default").head) else WdlString("")
          case Some(d) => throw new UnsupportedOperationException(s"Parameter ${v.variable} is required, but no value is specified")
          case None => throw new UnsupportedOperationException(s"This should not happen: could not find declaration for ${v.variable}")
        }
        case _ => throw new UnsupportedOperationException(s"Could not evaluate expression: ${expression.toString}")
      }
    }

    value match {
      case p: WdlPrimitive => p.valueString
      case a: WdlArray if attributes.contains("sep") => a.value.map(_.valueString).mkString(attributes.get("sep").head)
      case a: WdlArray => throw new UnsupportedOperationException(s"Expression '${expression.toString}' evaluated to an Array but no 'sep' was specified")
      case _ => throw new UnsupportedOperationException(s"Could not string-ify value: $value")
    }
  }
}
