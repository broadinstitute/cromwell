package cromwell.binding.command

import cromwell.binding.WdlExpression
import cromwell.binding.types.WdlType
import cromwell.binding.values.{WdlArray, WdlPrimitive, WdlString, WdlValue}

object ParameterCommandPart {
  val PostfixQuantifiersThatAcceptArrays = Set("+", "*")
  val OptionalPostfixQuantifiers = Set("?", "*")
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
    val paramValue = parameters.get(name) match {
      case Some(value) => value
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

    if (postfixQuantifier.isEmpty && wdlType != paramValueEvaluated.wdlType) {
      throw new UnsupportedOperationException(s"Incompatible type for $name: need a $wdlType, got a ${paramValue.wdlType}")
    }

    /* TODO: asString should be deprecated in the near future
     * It is being used as here because primitive types are trivially
     * turned into strings, but a more sophisticated solution will be
     * needed for compound types */
    paramValueEvaluated match {
      case param:WdlPrimitive => s"${prefix.getOrElse("")}${param.toWdlString}"
      case arr:WdlArray =>
        postfixQuantifier match {
          case Some(x) if ParameterCommandPart.PostfixQuantifiersThatAcceptArrays.contains(x) && attributes.contains("sep") =>
            val concatValue = arr.value.map {_.toWdlString}.mkString(attributes.get("sep").head)
            s"${prefix.getOrElse("")}$concatValue"
          case _ => throw new UnsupportedOperationException()
        }
      case _ => throw new UnsupportedOperationException("Not implemented yet")
    }
  }
}
