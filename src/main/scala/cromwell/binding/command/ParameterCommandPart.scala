package cromwell.binding.command

import cromwell.binding.WdlValue
import cromwell.binding.types.WdlType

case class ParameterCommandPart(wdlType: WdlType, name: String) extends CommandPart {
    override def toString: String = "${" + s"$wdlType $name" + "}"
    def instantiate(parameters: Map[String, WdlValue]): String = {
        val paramValue = parameters.getOrElse(name, {
            throw new UnsupportedOperationException(s"Parameter $name not found")
        })
        if (wdlType != paramValue.wdlType) {
            throw new UnsupportedOperationException(s"Incompatible type for $name: need a $wdlType, got a ${paramValue.wdlType}")
        }
        /* TODO: this also is a gross over-simplification */
        paramValue.value.toString
    }
}
