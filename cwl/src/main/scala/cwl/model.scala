package cwl

import cats.data.NonEmptyList
import eu.timepit.refined._
import cats.syntax.either._
import shapeless.{:+:, CNil, Witness}
import shapeless.syntax.singleton._
import cwl.LinkMergeMethod.LinkMergeMethod
import cwl.WorkflowStepInput.InputSource
import common.validation.ErrorOr.ErrorOr
import cwl.command.ParentName
import io.circe.Json
import wom.types.WomType
import wom.graph.WomIdentifier
import wom.graph.GraphNodePort.OutputPort
import wom.graph.expression.ExposedExpressionNode

case class WorkflowStepInput(
  id: String,
  source: Option[InputSource] = None,
  linkMerge: Option[LinkMergeMethod] = None,
  default: Option[CwlAny] = None,
  valueFrom: Option[StringOrExpression] = None) {

  def toExpressionNode(sourceMappings: Map[String, OutputPort],
                       outputTypeMap: Map[String, WomType],
                       inputs: Set[String]
                      )(implicit parentName: ParentName): ErrorOr[ExposedExpressionNode] = {
    val source = this.source.flatMap(_.select[String]).get
    val lookupId = FullyQualifiedName(source).id

    val outputTypeMapWithIDs = outputTypeMap.map {
      case (key, value) => FullyQualifiedName(key).id -> value
    }
    (for {
      inputType <- outputTypeMapWithIDs.get(lookupId).
        toRight(NonEmptyList.one(s"couldn't find $lookupId as derived from $source in map\n${outputTypeMapWithIDs.mkString("\n")}"))
      womExpression = WorkflowStepInputExpression(this, inputType, inputs)
      identifier = WomIdentifier(id)
      ret <- ExposedExpressionNode.fromInputMapping(identifier, womExpression, inputType, sourceMappings).toEither
    } yield ret).toValidated
  }
}

object WorkflowStepInput {
  type InputSource = String :+: Array[String] :+: CNil
}

object WorkflowStepInputSource {
  object String {
    def unapply(arg: InputSource): Option[String] = arg.select[String]
  }
  object StringArray {
    def unapply(arg: InputSource): Option[Array[String]] = arg.select[Array[String]]
  }
}

case class InputParameter(
                           id: String,
                           label: Option[String] = None,
                           secondaryFiles:
                             Option[
                               Expression :+:
                               String :+:
                               Array[
                                 Expression :+:
                                 String :+:
                                 CNil] :+:
                               CNil] = None,
                           format:
                             Option[
                               Expression :+:
                               String :+:
                               Array[String] :+:
                               CNil] = None,
                           streamable: Option[Boolean] = None,
                           doc: Option[String :+: Array[String] :+: CNil] = None,
                           inputBinding: Option[CommandLineBinding] = None,
                           default: Option[Json] = None, //can be of type "Any" which... sucks.
                           `type`: Option[MyriadInputType] = None) {

  type `type` = MyriadInputType
  type Id = String
}

case class InputRecordSchema(
  `type`: W.`"record"`.T,
  fields: Option[Array[InputRecordField]],
  label: Option[String])

case class InputRecordField(
  name: String,
  `type`: MyriadInputType,
  doc: Option[String],
  inputBinding: Option[CommandLineBinding],
  label: Option[String])

case class InputEnumSchema(
  symbols: Array[String],
  `type`: W.`"enum"`.T,
  label: Option[String],
  inputBinding: Option[CommandLineBinding])

case class InputArraySchema(
  items: MyriadInputType,
  `type`: W.`"array"`.T = Witness("array").value,
  label: Option[String] = None,
  inputBinding: Option[CommandLineBinding] = None)

case class CommandLineBinding(
                               loadContents: Option[Boolean] = None,
                               position: Option[Int] = None,
                               prefix: Option[String] = None,
                               separate: Option[String] = None,
                               itemSeparator: Option[String] = None,
                               valueFrom: Option[StringOrExpression] = None,
                               shellQuote: Option[Boolean] = None)

case class WorkflowOutputParameter(
                                    id: String,
                                    label: Option[String] = None,
                                    secondaryFiles:
                                      Option[
                                        Expression :+:
                                        String :+:
                                        Array[
                                          Expression :+:
                                          String :+:
                                          CNil] :+:
                                        CNil] = None,
                                    format: Option[Expression :+: String :+: Array[String] :+: CNil] = None,
                                    streamable: Option[Boolean] = None,
                                    doc: Option[String :+: Array[String] :+: CNil] = None,
                                    outputBinding: Option[CommandOutputBinding] = None,
                                    outputSource: Option[WorkflowOutputParameter#OutputSource] = None,
                                    linkMerge: Option[LinkMergeMethod] = None,
                                    `type`: Option[MyriadOutputType] = None) {

  type OutputSource = String :+: Array[String] :+: CNil
  type `type` = MyriadOutputType
  type Id = String
}

case class InputBinding(position: Int, prefix: String)

case class OutputRecordSchema(
  `type`: W.`"record"`.T,
  fields: Option[Array[OutputRecordField]],
  label: Option[String]
  )

case class OutputRecordField(
  name: String,
  `type`: MyriadOutputType,
  doc: Option[String],
  outputBinding: Option[CommandOutputBinding])

case class OutputEnumSchema(
  symbols: Array[String],
  `type`: W.`"enum"`.T,
  label: Option[String],
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
  `class`: W.`"SchemaDefRequirement"`.T,
  types: Array[InputRecordSchema :+: InputEnumSchema :+: InputArraySchema :+: CNil]
  )

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

case class InitialWorkDirRequirement(
  `class`: W.`"InitialWorkDirRequirement"`.T,
  listing:
    Array[
      File :+:
      Directory :+:
      Dirent :+:
      Expression :+:
      String :+:
      CNil
    ] :+:
    Expression :+:
    String :+:
    CNil)

/**
 *  Short for "Directory Entry"
 *  @see <a href="http://www.commonwl.org/v1.0/CommandLineTool.html#Dirent">Dirent Specification</a>
 */
case class Dirent(
                   entry: Expression :+: String :+: CNil,
                   entryName: Option[Expression :+: String :+: CNil],
                   writable: Option[Boolean])


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
                                coresMin: Long :+: Expression :+: String :+: CNil,
                                coresMax: Int :+: Expression :+: String :+: CNil,
                                ramMin: Long :+: Expression :+: String :+: CNil,
                                ramMax: Long :+: Expression :+: String :+: CNil,
                                tmpdirMin: Long :+: Expression :+: String :+: CNil,
                                tmpdirMax: Long :+: Expression :+: String :+: CNil,
                                outdirMin: Long :+: Expression :+: String :+: CNil,
                                outdirMax: Long :+: Expression :+: String :+: CNil)

case class SubworkflowFeatureRequirement(
  `class`: W.`"SubworkflowFeatureRequirement"`.T)

case class ScatterFeatureRequirement(
  `class`: W.`"ScatterFeatureRequirement"`.T)

case class MultipleInputFeatureRequirement(
  `class`: W.`"MultipleInputFeatureRequirement"`.T)

case class StepInputExpressionRequirement(
  `class`: W.`"StepInputExpressionRequirement"`.T)
