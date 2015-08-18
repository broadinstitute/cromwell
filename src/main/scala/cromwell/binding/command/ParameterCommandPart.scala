package cromwell.binding.command

import cromwell.binding.AstTools.EnhancedAstNode
import cromwell.binding.types.WdlType
import cromwell.binding.values.{WdlArray, WdlPrimitive, WdlString, WdlValue}
import cromwell.binding.{WdlExpression, WdlSyntaxErrorFormatter}
import cromwell.parser.WdlParser.{Ast, AstList, SyntaxError, Terminal}
import scala.language.postfixOps
import scala.collection.JavaConverters._

object ParameterCommandPart {
  val PostfixQuantifiersThatAcceptArrays = Set("+", "*")
  val OptionalPostfixQuantifiers = Set("?", "*")

  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): ParameterCommandPart = {
    val wdlType = ast.getAttribute("type").wdlType(wdlSyntaxErrorFormatter)
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    val prefix = Option(ast.getAttribute("prefix")) map { case t:Terminal => t.sourceString }
    val attributes = ast.getAttribute("attributes").asInstanceOf[AstList].asScala.toSeq map { a =>
      val ast = a.asInstanceOf[Ast]
      (ast.getAttribute("key").sourceString, ast.getAttribute("value").sourceString)
    } toMap
    val postfixQuantifier = ast.getAttribute("postfix") match {
      case t: Terminal => Some(t.sourceString)
      case _ => None
    }
    postfixQuantifier.foreach { quantifier =>
      val postfixQuantifierToken = ast.getAttribute("postfix").asInstanceOf[Terminal]
      if (PostfixQuantifiersThatAcceptArrays.contains(quantifier) && !attributes.contains("sep"))
        throw new SyntaxError(wdlSyntaxErrorFormatter.postfixQualifierRequiresSeparator(postfixQuantifierToken))
      if (!OptionalPostfixQuantifiers.contains(quantifier) && attributes.contains("default"))
        throw new SyntaxError(wdlSyntaxErrorFormatter.defaultAttributeOnlyAllowedForOptionalParameters(postfixQuantifierToken))
    }
    if (attributes.contains("default") && postfixQuantifier.isEmpty) {
      /* Calling .get below because we're assuming if attribute.contains("default"), then this operation is safe */
      val declarationOfDefaultAttribute = ast.getAttribute("attributes").asInstanceOf[AstList].asScala.toSeq.find{node =>
       node.asInstanceOf[Ast].getAttribute("key").asInstanceOf[Terminal].getSourceString == "default"
      }.get
      val terminal = declarationOfDefaultAttribute.asInstanceOf[Ast].getAttribute("key").asInstanceOf[Terminal]
      throw new SyntaxError(wdlSyntaxErrorFormatter.defaultAttributeOnlyAllowedForOptionalParameters(terminal))
    }
    new ParameterCommandPart(wdlType, name, prefix, attributes, postfixQuantifier)
  }
}

/**
 * Represents a parameter within a command, e.g. `${default="/etc/foo.conf" "--conf=" File conf?}`
 *
 * @param wdlType - The type that the input is required to be (`File`)
 * @param name - The name of the parameter (`conf`)
 * @param prefix - Optional prefix for when the parameter is specified.  This will be prepended to the parameter (`--conf=`)
 * @param attributes - A set of key (String) -> value (String) pairs from x="y" pairs within the parameter (default -> /etc/foo.conf)
 *                     There are only a subset of these attributes that have any meaning.  `default` is interpreted to be
 *                     the default value of the variable `conf` if no value is specified (since it's optional)
 * @param postfixQuantifier - The `?`, `*`, or `+` after the variable name.  This means "optional", "0-or-more", and "1-or-more", respectively
 */
case class ParameterCommandPart(wdlType: WdlType, name: String,
                                prefix: Option[String], attributes: Map[String, String],
                                postfixQuantifier: Option[String] = None) extends CommandPart {
  override def toString: String = "${" + s"${wdlType.toWdlString} $name" + "}"

  def instantiate(parameters: Map[String, WdlValue]): String = {
    def wasDefaultValueUsed: Boolean = parameters.get(name).isEmpty && attributes.contains("default")

    /* Order DOES matter in this match clause */
    val paramValue = parameters.get(name) match {
      case Some(value) => value
      case None if attributes.contains("default") => WdlString(attributes.get("default").get)
      case None if postfixQuantifier.isDefined && ParameterCommandPart.OptionalPostfixQuantifiers.contains(postfixQuantifier.head) => WdlString("")
      case _ => throw new UnsupportedOperationException(s"Parameter $name not found")
    }
    val paramValueEvaluated = paramValue match {
      case e: WdlExpression =>
        /* TODO: this should never happen because expressions will be evaluated
           before this method is called.  The reason why we can't always do it
           here is because unless scope information is stored in the WdlExpression
           itself (along with lookup and WdlFunctions objects), we'd have no way to
           evaluate this with any reasonable variable dereferencing.
         */
        throw new UnsupportedOperationException("All WdlExpressions must be evaluated before calling this")
      case x => x
    }

    if (!wasDefaultValueUsed && postfixQuantifier.isEmpty && wdlType != paramValueEvaluated.wdlType) {
      throw new UnsupportedOperationException(s"Incompatible type for $name: need a $wdlType, got a ${paramValue.wdlType}")
    }

    /* TODO: asString should be deprecated in the near future
     * It is being used as here because primitive types are trivially
     * turned into strings, but a more sophisticated solution will be
     * needed for compound types */
    paramValueEvaluated match {
      case WdlString("") => ""
      case param:WdlPrimitive => s"${prefix.getOrElse("")}${param.valueString}"
      case arr:WdlArray =>
        postfixQuantifier match {
          case Some(x) if ParameterCommandPart.PostfixQuantifiersThatAcceptArrays.contains(x) && attributes.contains("sep") =>
            val concatValue = arr.value.map {_.valueString}.mkString(attributes.get("sep").head)
            s"${prefix.getOrElse("")}$concatValue"
          case _ => throw new UnsupportedOperationException()
        }
      case _ => throw new UnsupportedOperationException("Not implemented yet")
    }
  }
}
