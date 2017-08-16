package wdl4s.cwl

import eu.timepit.refined._
import shapeless.{:+:, CNil, Witness}
import shapeless.syntax.singleton._
import wdl4s.cwl.CommandLineTool.StringOrExpression
import wdl4s.cwl.CommandOutputBinding.Glob
import wdl4s.cwl.LinkMergeMethod.LinkMergeMethod
import wdl4s.cwl.WorkflowStepInput.InputSource

case class WorkflowStepInput(
  id: String,
  source: Option[InputSource] = None,
  linkMerge: Option[LinkMergeMethod] = None,
  default: Option[CwlAny] = None,
  valueFrom: Option[ECMAScriptExpression :+: String :+: CNil] = None)

object WorkflowStepInput {
  type InputSource = String :+: Array[String] :+: CNil
}

case class InputParameter(
                           id: String,
                           label: Option[String] = None,
                           secondaryFiles:
                             Option[
                               ECMAScriptExpression :+:
                               String :+:
                               Array[
                                 ECMAScriptExpression :+:
                                 String :+:
                                 CNil] :+:
                               CNil] = None,
                           format:
                             Option[
                               ECMAScriptExpression :+:
                               String :+:
                               Array[String] :+:
                               CNil] = None,
                           streamable: Option[Boolean] = None,
                           doc: Option[String :+: Array[String] :+: CNil] = None,
                           inputBinding: Option[CommandLineBinding] = None,
                           default: Option[String] = None, //can be of type "Any" which... sucks.
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
  `type`: W.`"array"`.T,
  label: Option[String],
  inputBinding: Option[CommandLineBinding])

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
                                        ECMAScriptExpression :+:
                                        String :+:
                                        Array[
                                          ECMAScriptExpression :+:
                                          String :+:
                                          CNil] :+:
                                        CNil] = None,
                                    format: Option[ECMAScriptExpression :+: String :+: Array[String] :+: CNil] = None,
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

/** @see <a href="http://www.commonwl.org/v1.0/Workflow.html#CommandOutputBinding">CommandOutputBinding</a> */
case class CommandOutputBinding(
                                 glob: Option[Glob] = None,
                                 loadContents: Option[Boolean] = None,
                                 outputEval: Option[ECMAScriptExpression :+: String :+: CNil] = None)
object CommandOutputBinding {
  type Glob = ECMAScriptExpression :+: String :+: Array[String] :+: CNil

}

case class OutputArraySchema(
  items: MyriadOutputType,
  `type`: W.`"array"`.T,
  label: Option[String],
  outputBinding: Option[CommandOutputBinding])


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
      ECMAScriptExpression :+:
      String :+:
      CNil
    ] :+:
    ECMAScriptExpression :+:
    String :+:
    CNil)

/**
 *  Short for "Directory Entry"
 *  @see <a href="http://www.commonwl.org/v1.0/CommandLineTool.html#Dirent">Dirent Specification</a>
 */
case class Dirent(
                   entry: ECMAScriptExpression :+: String :+: CNil,
                   entryName: Option[ECMAScriptExpression :+: String :+: CNil],
                   writable: Option[Boolean])


case class EnvVarRequirement(
                              `class`: EnvVarRequirement.`class`.type = EnvVarRequirement.`class`,
                              envDef: Array[EnvironmentDef]
                            )

object EnvVarRequirement {
  val `class` : Witness.`"EnvVarRequirement"`.T = "EnvVarRequirement".narrow
}

case class EnvironmentDef(envName: String, envValue: StringOrExpression) {
  type EnvName = String
  type EnvValue = String
}


case class ShellCommandRequirement(`class`: W.`"ShellCommandRequirement"`.T = "ShellCommandRequirement".narrow)

case class ResourceRequirement(
                                `class`: W.`"ResourceRequirement"`.T,
                                coresMin: Long :+: ECMAScriptExpression :+: String :+: CNil,
                                coresMax: Int :+: ECMAScriptExpression :+: String :+: CNil,
                                ramMin: Long :+: ECMAScriptExpression :+: String :+: CNil,
                                ramMax: Long :+: ECMAScriptExpression :+: String :+: CNil,
                                tmpdirMin: Long :+: ECMAScriptExpression :+: String :+: CNil,
                                tmpdirMax: Long :+: ECMAScriptExpression :+: String :+: CNil,
                                outdirMin: Long :+: ECMAScriptExpression :+: String :+: CNil,
                                outdirMax: Long :+: ECMAScriptExpression :+: String :+: CNil)

case class SubworkflowFeatureRequirement(
  `class`: W.`"SubworkflowFeatureRequirement"`.T)

case class ScatterFeatureRequirement(
  `class`: W.`"ScatterFeatureRequirement"`.T)

case class MultipleInputFeatureRequirement(
  `class`: W.`"MultipleInputFeatureRequirement"`.T)

case class StepInputExpressionRequirement(
  `class`: W.`"StepInputExpressionRequirement"`.T)
