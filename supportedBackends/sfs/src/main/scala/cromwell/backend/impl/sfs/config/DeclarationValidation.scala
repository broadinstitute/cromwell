package cromwell.backend.impl.sfs.config

import cromwell.backend.validation._
import wdl4s.expression.NoFunctions
import wdl4s.types._
import wdl4s.values.WdlValue
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
      case name if name == DockerValidation.instance.key =>
        new DeclarationValidation(declaration, DockerValidation.instance)
      case name if name == CpuValidation.instance.key => new DeclarationValidation(declaration, CpuValidation.instance)
      // See MemoryDeclarationValidation for more info
      case name if MemoryDeclarationValidation.isMemoryDeclaration(name) =>
        new MemoryDeclarationValidation(declaration)
      // All other declarations must be a Boolean, Float, Integer, or String.
      case _ =>
        val validatedRuntimeAttr = validator(declaration.wdlType, declaration.unqualifiedName)
        new DeclarationValidation(declaration, validatedRuntimeAttr)
    }
  }

  private def validator(wdlType: WdlType, unqualifiedName: String): PrimitiveRuntimeAttributesValidation[_, _] = {
    wdlType match {
      case WdlBooleanType => new BooleanRuntimeAttributesValidation(unqualifiedName)
      case WdlFloatType => new FloatRuntimeAttributesValidation(unqualifiedName)
      case WdlIntegerType => new IntRuntimeAttributesValidation(unqualifiedName)
      case WdlStringType => new StringRuntimeAttributesValidation(unqualifiedName)
      case WdlOptionalType(x) => validator(x, unqualifiedName)
      case other => throw new RuntimeException(s"Unsupported config runtime attribute $other $unqualifiedName")
    }
  }
}

/**
  * Tracks some declaration with a validation for the declaration. Starts with a basic instance of a validation.
  *
  * @param declaration        The declaration from the config "runtime-attributes".
  * @param instanceValidation A basic instance validation for the declaration.
  */
class DeclarationValidation(declaration: Declaration, instanceValidation: RuntimeAttributesValidation[_]) {
  val key: String = declaration.unqualifiedName

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
    validation.withDefault(Some(wdlExpression.evaluate(NoLookup, NoFunctions).get))
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
