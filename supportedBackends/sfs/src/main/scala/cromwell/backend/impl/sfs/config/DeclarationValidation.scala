package cromwell.backend.impl.sfs.config

import cromwell.backend.impl.sfs.config.ConfigConstants._
import cromwell.backend.validation._
import wdl.draft2.model.expression.NoFunctions
import wdl.draft2.model.{Declaration, NoLookup, WdlExpression}
import wom.RuntimeAttributesKeys
import wom.types._
import wom.values.WomValue

import scala.annotation.tailrec

/**
  * Creates instances of runtime attribute validations from WDL declarations.
  */
object DeclarationValidation {
  def fromDeclarations(declarations: Seq[Declaration], callCachedRuntimeAttributes: Map[String, Boolean]): Seq[DeclarationValidation] = {
    declarations map fromDeclaration(callCachedRuntimeAttributesMap = callCachedRuntimeAttributes) _
  }

  /**
    * Create a runtime attribute validation from a WDL declaration.
    *
    * @param declaration The declaration.
    * @return The DeclarationValidation object for the declaration.
    */
  def fromDeclaration(callCachedRuntimeAttributesMap: Map[String, Boolean])(declaration: Declaration): DeclarationValidation = {
    declaration.unqualifiedName match {
      // Docker and CPU are special keys understood by cromwell.
      case name if name == DockerValidation.instance.key =>
        new DeclarationValidation(declaration, DockerValidation.instance, usedInCallCachingOverride = None)
      case RuntimeAttributesKeys.CpuKey => new CpuDeclarationValidation(declaration, CpuValidation.instance)
      case RuntimeAttributesKeys.CpuMinKey => new CpuDeclarationValidation(declaration, CpuValidation.instanceMin)
      case RuntimeAttributesKeys.CpuMaxKey => new CpuDeclarationValidation(declaration, CpuValidation.instanceMax)
      // See MemoryDeclarationValidation for more info
      case name if MemoryDeclarationValidation.isMemoryDeclaration(name, MemoryRuntimeAttribute, MemoryRuntimeAttributePrefix) =>
        new MemoryDeclarationValidation(declaration, MemoryRuntimeAttribute, MemoryRuntimeAttributePrefix)
      case name if MemoryDeclarationValidation.isMemoryDeclaration(name, MemoryMinRuntimeAttribute, MemoryRuntimeAttributePrefix) =>
        new MemoryDeclarationValidation(declaration, MemoryMinRuntimeAttribute, MemoryMinRuntimeAttributePrefix)
      case name if MemoryDeclarationValidation.isMemoryDeclaration(name, MemoryMaxRuntimeAttribute, MemoryRuntimeAttributePrefix) =>
        new MemoryDeclarationValidation(declaration, MemoryMaxRuntimeAttribute, MemoryMaxRuntimeAttributePrefix)
      case name if MemoryDeclarationValidation.isMemoryDeclaration(name, DiskRuntimeAttribute, DiskRuntimeAttributePrefix) =>
        new MemoryDeclarationValidation(declaration, DiskRuntimeAttribute, DiskRuntimeAttributePrefix)
      // All other declarations must be a Boolean, Float, Integer, or String.
      case _ =>
        val validatedRuntimeAttr = validator(declaration.womType, declaration.unqualifiedName)
        new DeclarationValidation(
          declaration = declaration,
          instanceValidation = validatedRuntimeAttr,
          usedInCallCachingOverride = callCachedRuntimeAttributesMap.get(declaration.unqualifiedName)
        )
    }
  }

  @tailrec
  private def validator(womType: WomType, unqualifiedName: String): PrimitiveRuntimeAttributesValidation[_, _] = {
    womType match {
      case WomBooleanType => new BooleanRuntimeAttributesValidation(unqualifiedName)
      case WomFloatType => new FloatRuntimeAttributesValidation(unqualifiedName)
      case WomIntegerType => new IntRuntimeAttributesValidation(unqualifiedName)
      case WomStringType => new StringRuntimeAttributesValidation(unqualifiedName)
      case WomOptionalType(x) => validator(x, unqualifiedName)
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
class DeclarationValidation(declaration: Declaration, instanceValidation: RuntimeAttributesValidation[_], usedInCallCachingOverride: Option[Boolean]) {
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
    import scala.language.existentials

    // The RuntimeAttributesValidation object contains functions to wrap existing validation functions into new validations with modified behaviors.
    // As a first approximation, think "caseClass.copy, but for validation functions"
    // In this case, we might (or might not) want to make our validations:
    // 1. have defaults:
    val validationWithDefault = if (declaration.expression.isDefined) default(instanceValidation, declaration.expression.get) else instanceValidation
    // 2. be optional:
    val validationWithDefaultAndOptionality = if (declaration.womType.isInstanceOf[WomOptionalType]) validationWithDefault.optional else validationWithDefault
    // Or 3. have customized call caching properties:
    val validationWithDefaultAndOptionalityAndCallCaching = usedInCallCachingOverride match {
      case Some(usedInCallCachingValue) => RuntimeAttributesValidation.withUsedInCallCaching(validationWithDefaultAndOptionality, usedInCallCachingValue)
      case None => validationWithDefaultAndOptionality
    }

    validationWithDefaultAndOptionalityAndCallCaching
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
  def extractWdlValueOption(validatedRuntimeAttributes: ValidatedRuntimeAttributes): Option[WomValue] = {
    RuntimeAttributesValidation.extractOption(instanceValidation, validatedRuntimeAttributes) map {
      declaration.womType.coerceRawValue(_).get
    }
  }
}
