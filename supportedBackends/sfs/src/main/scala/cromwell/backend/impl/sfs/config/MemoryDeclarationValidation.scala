package cromwell.backend.impl.sfs.config

import cromwell.backend.MemorySize
import cromwell.backend.validation._
import wdl4s.wdl.expression.NoFunctions
import wdl4s.parser.MemoryUnit
import wdl4s.wdl.types._
import wdl4s.wdl.values.{WdlFloat, WdlInteger, WdlOptionalValue, WdlValue}
import wdl4s.wdl.{Declaration, NoLookup, WdlExpression}

/**
  * Maps declarations of memory in WDL runtime attributes to the commands used to submit.
  *
  * The wdl runtime attributes for memory specified as strings such as:
  *
  * {{{
  *   runtime {
  *     memory: "500 MB"
  *   }
  * }}}
  *
  * However, the backend configuration only supports specifying memory in amounts in Float or Int. To specify the unit
  * of the amount, the string "memory_" is suffixed with the unit, for example "memory_mb" or "memory_gb", or even
  * "memory_ki".
  *
  * This class and companion object will do the conversion. The backend configuration should use the runtime attribute
  * "Float? memory_gb", meaning that "memory" is now an optional runtime attribute, and will be converted to GB.
  *
  * Just like the runtime attribute, when no units are specified the config, the default unit is bytes.
  *
  * @param declaration The declaration used to create this memory validation.
  */
class MemoryDeclarationValidation(declaration: Declaration, attributeName: String, attributeNamePrefix: String)
  extends DeclarationValidation(declaration, MemoryValidation.instance(attributeName)) {

  import MemoryDeclarationValidation._

  /**
    * Converts the validation to a version with the default from the memory expression.
    *
    * If the backend configuration contains a runtime attribute such as "Float memory_gb = 1.0", then the default will
    * be set to 1 GB of memory when the attribute is not set.
    *
    * @param validation    The validation to set the default for.
    * @param wdlExpression The declaration expression to retrieve the default.
    * @return The new validation.
    */
  override protected def default(validation: RuntimeAttributesValidation[_],
                                 wdlExpression: WdlExpression): RuntimeAttributesValidation[_] = {
    val wdlValue = declaration.expression.get.evaluate(NoLookup, NoFunctions).get
    val amount: Double = defaultAmount(wdlValue)
    val memorySize = MemorySize(amount, declarationMemoryUnit)
    validation.withDefault(WdlInteger(memorySize.bytes.toInt))
  }

  private def defaultAmount(wdlValue: WdlValue): Double = {
    wdlValue match {
      case WdlInteger(value) => value.toDouble
      case WdlFloat(value) => value
      case WdlOptionalValue(_, Some(optionalWdlValue)) => defaultAmount(optionalWdlValue)
      case other => throw new RuntimeException(s"Unsupported memory default: $other")
    }
  }

  private lazy val declarationMemoryUnit: MemoryUnit = {
    val suffix = memoryUnitSuffix(declaration.unqualifiedName, attributeName, attributeNamePrefix)
    val memoryUnitOption = MemoryUnit.values.find(_.suffixes.map(_.toLowerCase).contains(suffix.toLowerCase))
    memoryUnitOption match {
      case Some(memoryUnit) => memoryUnit
      case None => throw new IllegalArgumentException(s"MemoryUnit with suffix $suffix was not found.")
    }
  }

  /**
    * Converts the memory value from a `MemorySize` to a `Float` or `Int` based on the units.
    *
    * @param validatedRuntimeAttributes The validated attributes.
    * @return The value from the collection wrapped in `Some`, or `None` if the value wasn't found.
    */
  override def extractWdlValueOption(validatedRuntimeAttributes: ValidatedRuntimeAttributes): Option[WdlValue] = {
    RuntimeAttributesValidation.extractOption(MemoryValidation.instance(attributeName), validatedRuntimeAttributes) map
      coerceMemorySize(declaration.wdlType)
  }

  private def coerceMemorySize(wdlType: WdlType)(value: MemorySize): WdlValue = {
    wdlType match {
      case WdlIntegerType => WdlInteger(value.to(declarationMemoryUnit).amount.toInt)
      case WdlFloatType => WdlFloat(value.to(declarationMemoryUnit).amount)
      case WdlOptionalType(optionalType) => coerceMemorySize(optionalType)(value)
      case other => throw new RuntimeException(s"Unsupported wdl type for memory: $other")
    }
  }
}

object MemoryDeclarationValidation {
  def isMemoryDeclaration(name: String, attributeName: String, attributeNamePrefix: String): Boolean = {
    name match {
      case `attributeName` => true
      case prefixed if prefixed.startsWith(attributeNamePrefix) =>
        val suffix = memoryUnitSuffix(name, attributeName, attributeNamePrefix)
        MemoryUnit.values exists {
          _.suffixes.map(_.toLowerCase).contains(suffix)
        }
      case _ => false
    }
  }

  private def memoryUnitSuffix(name: String, attributeName: String, attributeNamePrefix: String) = {
    if (name == attributeName)
      MemoryUnit.Bytes.suffixes.head
    else
      name.substring(attributeNamePrefix.length)
  }
}
