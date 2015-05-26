package cromwell.binding.command

import cromwell.binding.types.WdlType
import cromwell.binding.values.{WdlPrimitive, WdlValue}

case class ParameterCommandPart(wdlType: WdlType, name: String) extends CommandPart {
  override def toString: String = "${" + s"$wdlType $name" + "}"

  def instantiate(parameters: Map[String, WdlValue]): String = {
    val paramValue = parameters.getOrElse(name, {
      throw new UnsupportedOperationException(s"Parameter $name not found")
    })
    if (wdlType != paramValue.wdlType) {
      throw new UnsupportedOperationException(s"Incompatible type for $name: need a $wdlType, got a ${paramValue.wdlType}")
    }

    /* TODO: asString should be deprecated in the near future
     * It is being used as here because primitive types are trivially
     * turned into strings, but a more sophisticated solution will be
     * needed for compound types */
    paramValue match {
      case param:WdlPrimitive => param.asString
      case _ => throw new UnsupportedOperationException("Not implemented yet")
    }
  }
}
