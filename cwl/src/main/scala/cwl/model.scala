package cwl

import cwl.CommandLineTool.{CommandBindingSortingKey, SortKeyAndCommandPart}
import cwl.SchemaDefRequirement.{AsInputEnumSchema, AsInputRecordSchema, SchemaDefTypes}
import cwl.CwlType.CwlType
import cwl.WorkflowStepInput.InputSource
import cwl.command.ParentName
import cwl.internal.GigabytesToBytes
import eu.timepit.refined._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.string.MatchesRegex
import shapeless.syntax.singleton._
import shapeless.{:+:, CNil, Coproduct, Poly1, Witness}
import wom.types.WomType
import wom.values.WomValue
import mouse.all._

object WorkflowStepInputSource {
  object String {
    def unapply(arg: InputSource): Option[String] = arg.select[String]
  }
  object StringArray {
    def unapply(arg: InputSource): Option[Array[String]] = arg.select[Array[String]]
  }
}

/**
  * Describes a bespoke type.
  *
  * @param name This field actually does _not_ appear in the v1.0 schema, but it is used anyway in the conformance tests.
  *             After some consideration it was determined that we should close our eyes and pretend it is in the spec.  It
  *             makes its formal appearance as a required field in v1.1.
  */
case class InputRecordSchema(
  name: String,
  fields: Option[Array[InputRecordField]] = None,
  `type`: W.`"record"`.T = W("record").value,
  label: Option[String] = None) {

}

case class InputRecordField(
  name: String,
  `type`: MyriadInputType,
  doc: Option[String] = None,
  inputBinding: Option[InputCommandLineBinding],
  label: Option[String] = None)

case class InputArraySchema
(
  items: MyriadInputType,
  `type`: W.`"array"`.T = Witness("array").value,
  label: Option[String] = None,
  inputBinding: Option[InputCommandLineBinding] = None,
  // IAS.secondaryFiles are NOT listed in 1.0 spec, but according to jgentry they will be, maybe
  secondaryFiles: Option[SecondaryFiles] = None
)

trait CommandLineBinding {
  def loadContents: Option[Boolean]
  def position: Option[Int]
  def prefix: Option[String]
  def separate: Option[Boolean]
  def itemSeparator: Option[String]
  def optionalValueFrom: Option[StringOrExpression]
  def shellQuote: Option[Boolean]
  // separate defaults to true
  def effectiveSeparate = separate.getOrElse(true)
}

object InputCommandLineBinding {
  def default = InputCommandLineBinding()
}

case class InputCommandLineBinding(
                               loadContents: Option[Boolean] = None,
                               position: Option[Int] = None,
                               prefix: Option[String] = None,
                               separate: Option[Boolean] = None,
                               itemSeparator: Option[String] = None,
                               valueFrom: Option[StringOrExpression] = None,
                               shellQuote: Option[Boolean] = None) extends CommandLineBinding {
  override val optionalValueFrom = valueFrom

  def toCommandPart(sortingKey: CommandBindingSortingKey, boundValue: WomValue, hasShellCommandRequirement: Boolean, expressionLib: ExpressionLib) = {
    SortKeyAndCommandPart(sortingKey, InputCommandLineBindingCommandPart(this, boundValue)(hasShellCommandRequirement, expressionLib))
  }
}

// valueFrom is required for command line bindings in the argument section: http://www.commonwl.org/v1.0/CommandLineTool.html#CommandLineBinding
case class ArgumentCommandLineBinding(
                               valueFrom: StringOrExpression,
                               loadContents: Option[Boolean] = None,
                               position: Option[Int] = None,
                               prefix: Option[String] = None,
                               separate: Option[Boolean] = None,
                               itemSeparator: Option[String] = None,
                               shellQuote: Option[Boolean] = None) extends CommandLineBinding {
  override val optionalValueFrom = Option(valueFrom)
}

case class InputBinding(position: Int, prefix: String)

object MyriadOutputInnerTypeCacheableString extends Poly1 {
  import Case._

  private def cacheableOutputRecordFieldString(field: OutputRecordField): String = {
    val fieldType = field.`type`.fold(MyriadOutputTypeCacheableString)
    val lqn = field.name.substring(field.name.lastIndexOf('#') + 1)
    s"OutputRecordField($lqn,$fieldType,${field.doc},${field.outputBinding})"
  }

  implicit def recordSchema: Aux[OutputRecordSchema, String] = at[OutputRecordSchema] {
    s =>
      val t = s.`type`
      s"OutputRecordSchema($t,${s.fields.map(a => "Array(" + a.map(cacheableOutputRecordFieldString).mkString(",") + ")" )},${s.label})"
  }

  implicit def arraySchema: Aux[OutputArraySchema, String] = at[OutputArraySchema] {
    a =>
      val is: String = a.items.fold(MyriadOutputTypeCacheableString)
      val t = a.`type`
      s"OutputArraySchema($is,$t,${a.label},${a.outputBinding})"
  }

  implicit def enumSchema: Aux[OutputEnumSchema, String] = at[OutputEnumSchema] { _.toString }
  implicit def cwlType: Aux[CwlType, String] = at[CwlType] { _.toString }
  implicit def string: Aux[String, String] = at[String] { identity }
}


object MyriadOutputTypeCacheableString extends Poly1 {
  import Case._

  implicit def one: Aux[MyriadOutputInnerType, String] = at[MyriadOutputInnerType] {
    _.fold(MyriadOutputInnerTypeCacheableString)
  }

  implicit def many: Aux[Array[MyriadOutputInnerType], String] = at[Array[MyriadOutputInnerType]] { a =>
    val strings: Array[String] = a.map(_.fold(MyriadOutputInnerTypeCacheableString))
    "Array(" + strings.mkString(",") + ")"
  }
}

object SecondaryFilesCacheableString extends Poly1 {
  import Case._

  implicit def one: Aux[StringOrExpression, String] = at[StringOrExpression] {
    _.toString
  }

  implicit def array: Aux[Array[StringOrExpression], String] = at[Array[StringOrExpression]] {
    _.mkString("Array(", ",", ")")
  }

}

case class OutputRecordSchema(
  `type`: W.`"record"`.T,
  fields: Option[Array[OutputRecordField]],
  label: Option[String])

case class OutputRecordField(
  name: String,
  `type`: MyriadOutputType,
  doc: Option[String],
  outputBinding: Option[CommandOutputBinding])

case class OutputArraySchema(
  items: MyriadOutputType,
  `type`: W.`"array"`.T = Witness("array").value,
  label: Option[String] = None,
  outputBinding: Option[CommandOutputBinding] = None)

case class InlineJavascriptRequirement(
  `class`: W.`"InlineJavascriptRequirement"`.T = "InlineJavascriptRequirement".narrow,
  expressionLib: Option[Array[String]] = None)

case class SchemaDefRequirement(
  types: Array[SchemaDefTypes] = Array.empty,
  `class`: W.`"SchemaDefRequirement"`.T = Witness("SchemaDefRequirement").value) {

  def lookupType(tpe: String): Option[WomType] =
    lookupCwlType(tpe).flatMap{
      case AsInputRecordSchema(inputRecordSchema: InputRecordSchema) => MyriadInputInnerTypeToWomType.inputRecordSchemaToWomType(inputRecordSchema).apply(this) |> Option.apply
      case _ => None
    }

  //Currently only InputRecordSchema has a name in the spec, so it is the only thing that can be referenced via string
  def lookupCwlType(tpe: String): Option[SchemaDefTypes] = {

    def matchesType(inputEnumSchema: InputEnumSchema): Boolean = {
      inputEnumSchema.name.fold(false){name => FileAndId(name)(ParentName.empty).id equalsIgnoreCase FileAndId(tpe)(ParentName.empty).id}
    }

    types.toList.flatMap {
      case AsInputRecordSchema(inputRecordSchema: InputRecordSchema) if FileAndId(inputRecordSchema.name)(ParentName.empty).id equalsIgnoreCase FileAndId(tpe)(ParentName.empty).id =>
        List(Coproduct[SchemaDefTypes](inputRecordSchema))
      case AsInputEnumSchema(inputEnumSchema: InputEnumSchema) if matchesType(inputEnumSchema) =>
        List(Coproduct[SchemaDefTypes](inputEnumSchema))
      case _ =>
        List()
    }.headOption
  }
}

object SchemaDefRequirement {
  type SchemaDefTypes = InputRecordSchema :+: InputEnumSchema :+: InputArraySchema :+: CNil

  object AsInputRecordSchema {
    def unapply(arg: SchemaDefTypes): Option[InputRecordSchema] = arg.select[InputRecordSchema]
  }

  object AsInputEnumSchema {
    def unapply(arg: SchemaDefTypes): Option[InputEnumSchema] = arg.select[InputEnumSchema]
  }

  object AsInputArraySchema {
    def unapply(arg: SchemaDefTypes): Option[InputArraySchema] = arg.select[InputArraySchema]
  }
}

//There is a large potential for regex refinements on these string types
case class DockerRequirement(
  `class`: W.`"DockerRequirement"`.T,
  dockerPull: Option[String], //TODO Refine to match a docker image regex?
  dockerLoad: Option[String],
  dockerFile: Option[String],
  dockerImport: Option[String],
  dockerImageId: Option[String],
  dockerOutputDirectory: Option[String]
  )

case class SoftwareRequirement(
  `class`: W.`"SoftwareRequirement"`.T,
  packages: Array[SoftwarePackage] = Array.empty
  )

case class SoftwarePackage(
  `package`: String,
  version: Option[Array[String]],
  specs: Option[Array[String]] // This could be refined to match a regex for IRI.
  ) {
  type Package = String
  type Specs = Array[String]
}

case class EnvVarRequirement(
                              `class`: EnvVarRequirement.ClassType = EnvVarRequirement.`class`,
                              envDef: Array[EnvironmentDef]
                            )

object EnvVarRequirement {
  type ClassType = Witness.`"EnvVarRequirement"`.T

  val `class`: ClassType = "EnvVarRequirement".asInstanceOf[ClassType]
}

case class EnvironmentDef(envName: String, envValue: StringOrExpression) {
  type EnvName = String
  type EnvValue = String
}


case class ShellCommandRequirement(`class`: W.`"ShellCommandRequirement"`.T = "ShellCommandRequirement".narrow)

case class ResourceRequirement(
                                `class`: W.`"ResourceRequirement"`.T,
                                coresMin: Option[ResourceRequirementType],
                                coresMax: Option[ResourceRequirementType],
                                ramMin: Option[ResourceRequirementType],
                                ramMax: Option[ResourceRequirementType],
                                tmpdirMin: Option[ResourceRequirementType],
                                tmpdirMax: Option[ResourceRequirementType],
                                outdirMin: Option[ResourceRequirementType],
                                outdirMax: Option[ResourceRequirementType]) {
  def effectiveCoreMin = coresMin.orElse(coresMax)
  def effectiveCoreMax = coresMax.orElse(coresMin)

  def effectiveRamMin = ramMin.orElse(ramMax).map(_.fold(GigabytesToBytes))
  def effectiveRamMax = ramMax.orElse(ramMin).map(_.fold(GigabytesToBytes))

  def effectiveTmpdirMin = tmpdirMin.orElse(tmpdirMax)
  def effectiveTmpdirMax = tmpdirMax.orElse(tmpdirMin)

  def effectiveOutdirMin = outdirMin.orElse(outdirMax)
  def effectiveOutdirMax = outdirMax.orElse(outdirMin)
}

/**
  * This promotes DNA Nexus InputResourceRequirement to a first class citizen requirement, which it really isn't.
  * Since it's the only one for now it's not a big deal but if more of these pop up we might want to treat custom requirements
  * in a different way
  */
case class DnaNexusInputResourceRequirement(
                                     `class`: String Refined MatchesRegex[W.`".*InputResourceRequirement"`.T],
                                     indirMin: Option[Long]
                                   )

case class SubworkflowFeatureRequirement(
  `class`: W.`"SubworkflowFeatureRequirement"`.T)

case class ScatterFeatureRequirement(
  `class`: W.`"ScatterFeatureRequirement"`.T)

case class MultipleInputFeatureRequirement(
  `class`: W.`"MultipleInputFeatureRequirement"`.T)

case class StepInputExpressionRequirement(
  `class`: W.`"StepInputExpressionRequirement"`.T)
