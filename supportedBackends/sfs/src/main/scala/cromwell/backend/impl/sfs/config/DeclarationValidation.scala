package cromwell.backend.impl.sfs.config

import cromwell.backend.MemorySize
import cromwell.backend.impl.sfs.config.ConfigConstants._
import cromwell.backend.validation._
import wdl4s.expression.NoFunctions
import wdl4s.parser.MemoryUnit
import wdl4s.types._
import wdl4s.values.{WdlFloat, WdlInteger, WdlValue}
import wdl4s.{Declaration, NoLookup, WdlExpression}

/**
  * Creates instances of runtime attribute validations from WDL declarations.
  */
object DeclarationValidation {
  def fromDeclarations(declarations: Seq[Declaration]): Seq[DeclarationValidation] = {
    declarations map fromDeclaration
  }

  /**
    * Create a runtime attribute validation from a WDL declaration.
    *
    * @param declaration The declaration.
    * @return The DeclarationValidation object for the declaration.
    */
  def fromDeclaration(declaration: Declaration): DeclarationValidation = {
    declaration.unqualifiedName match {
      // Docker and CPU are special keys understood by cromwell.
      case DockerValidation.key => new DeclarationValidation(declaration, DockerValidation.instance)
      case CpuValidation.key => new DeclarationValidation(declaration, CpuValidation.default)
      // See MemoryDeclarationValidation for more info
      case name if MemoryDeclarationValidation.isMemoryDeclaration(name) =>
        new MemoryDeclarationValidation(declaration)
      // All other declarations must be a Boolean, Float, Integer, or String.
      case _ =>
        val validatedRuntimeAttr = validator(declaration.wdlType, declaration.unqualifiedName)
        new DeclarationValidation(declaration, validatedRuntimeAttr)
    }
  }

  private def validator(wdlType: WdlType, unqualifiedName: String): PrimitiveRuntimeAttributesValidation[_] = wdlType match {
    case WdlBooleanType => new BooleanRuntimeAttributesValidation(unqualifiedName)
    case WdlFloatType => new FloatRuntimeAttributesValidation(unqualifiedName)
    case WdlIntegerType => new IntRuntimeAttributesValidation(unqualifiedName)
    case WdlStringType => new StringRuntimeAttributesValidation(unqualifiedName)
    case WdlOptionalType(x) => validator(x, unqualifiedName)
    case other => throw new RuntimeException(s"Unsupported config runtime attribute $other $unqualifiedName")
  }
}

/**
  * Tracks some declaration with a validation for the declaration. Starts with a basic instance of a validation.
  *
  * @param declaration        The declaration from the config "runtime-attributes".
  * @param instanceValidation A basic instance validation for the declaration.
  */
class DeclarationValidation(declaration: Declaration, instanceValidation: RuntimeAttributesValidation[_]) {
  val key = declaration.unqualifiedName

  /**
    * Creates a validation, by adding on defaults if they're specified in the declaration, and then making the
    * declaration optional if a "?" was placed on the declaration.
    *
    * Examples:
    * {{{
    *   String  custom_attribute             # Required, no default
    *   String? custom_attribute             # Optional, no default
    *   String  custom_attribute = "default" # Required, but defaulted to the value "default"
    *   String? custom_attribute = "default" # Optional, but defaulted to the value "default"
    * }}}
    *
    * The last one really shouldn't be used, but the syntax is technically valid, currently.
    *
    * @return The validation.
    */
  def makeValidation(): RuntimeAttributesValidation[_] = {
    val validationDefault = if (declaration.expression.isDefined)
      default(instanceValidation, declaration.expression.get)
    else instanceValidation
    if (declaration.wdlType.isInstanceOf[WdlOptionalType]) validationDefault.optional else validationDefault
  }

  /**
    * Overrides the passed in validation with a default from the wdl expression.
    *
    * @param validation    The original validation.
    * @param wdlExpression The wdl expression with the default value.
    * @return A new copy of the validation with the default value.
    */
  protected def default(validation: RuntimeAttributesValidation[_],
                        wdlExpression: WdlExpression): RuntimeAttributesValidation[_] = {
    validation.withDefault(wdlExpression.evaluate(NoLookup, NoFunctions).get)
  }

  /**
    * Utility to get the value of this declaration from a collection of validated runtime attributes.
    *
    * @param validatedRuntimeAttributes The validated attributes.
    * @return The value from the collection wrapped in `Some`, or `None` if the value wasn't found.
    */
  def extractWdlValueOption(validatedRuntimeAttributes: ValidatedRuntimeAttributes): Option[WdlValue] = {
    RuntimeAttributesValidation.extractOption(instanceValidation, validatedRuntimeAttributes) map {
      declaration.wdlType.coerceRawValue(_).get
    }
  }
}

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
class MemoryDeclarationValidation(declaration: Declaration)
  extends DeclarationValidation(declaration, MemoryValidation.instance) {

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
  override protected def default(validation: RuntimeAttributesValidation[_], wdlExpression: WdlExpression) = {
    val wdlValue = declaration.expression.get.evaluate(NoLookup, NoFunctions).get
    val amount: Double = wdlValue match {
      case WdlInteger(value) => value.toDouble
      case WdlFloat(value) => value
      case other => throw new RuntimeException(s"Unsupported memory default: $other")
    }
    val memorySize = MemorySize(amount, declarationMemoryUnit)
    validation.withDefault(WdlInteger(memorySize.bytes.toInt))
  }

  private lazy val declarationMemoryUnit: MemoryUnit = {
    val suffix = memoryUnitSuffix(declaration.unqualifiedName)
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
    RuntimeAttributesValidation.extractOption(MemoryValidation.instance, validatedRuntimeAttributes) map { value =>
      declaration.wdlType match {
        case WdlIntegerType => WdlInteger(value.to(declarationMemoryUnit).amount.toInt)
        case WdlFloatType => WdlFloat(value.to(declarationMemoryUnit).amount)
        case other => throw new RuntimeException(s"Unsupported wdl type for memory: $other")
      }
    }
  }
}

object MemoryDeclarationValidation {
  def isMemoryDeclaration(name: String): Boolean = {
    name match {
      case MemoryRuntimeAttribute => true
      case prefixed if prefixed.startsWith(MemoryRuntimeAttributePrefix) =>
        val suffix = memoryUnitSuffix(name)
        MemoryUnit.values exists {
          _.suffixes.map(_.toLowerCase).contains(suffix)
        }
      case _ => false
    }
  }

  private def memoryUnitSuffix(name: String) = {
    if (name == MemoryRuntimeAttribute)
      MemoryUnit.Bytes.suffixes.head
    else
      name.substring(MemoryRuntimeAttributePrefix.length)
  }
}
