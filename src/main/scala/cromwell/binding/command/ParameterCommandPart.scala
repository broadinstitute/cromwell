package cromwell.binding.command

import cromwell.binding.WdlExpression
import cromwell.binding.types.WdlType
import cromwell.binding.values.{WdlString, WdlArray, WdlPrimitive, WdlValue}
import scala.util.Success

case class ParameterCommandPart(wdlType: WdlType, name: String,
                                prefix: Option[String], attributes: Map[String, String],
                                postfixQuantifier: Option[String]) extends CommandPart {
  override def toString: String = "${" + s"$wdlType $name" + "}"

  def instantiate(parameters: Map[String, WdlValue]): String = {
    val paramValue = parameters.get(name) match {
      case Some(value) => value
      case None if postfixQuantifier.isDefined && Set("?", "*").contains(postfixQuantifier.head) => WdlString("")
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

    if (postfixQuantifier.isEmpty && wdlType != paramValueEvaluated.wdlType) {
      throw new UnsupportedOperationException(s"Incompatible type for $name: need a $wdlType, got a ${paramValue.wdlType}")
    }

    /* TODO: asString should be deprecated in the near future
     * It is being used as here because primitive types are trivially
     * turned into strings, but a more sophisticated solution will be
     * needed for compound types */
    paramValueEvaluated match {
      case param:WdlPrimitive => s"${prefix.getOrElse("")}${param.asString}"
      case arr:WdlArray => ""
      case _ => throw new UnsupportedOperationException("Not implemented yet")
    }
  }
}
