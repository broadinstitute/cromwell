package wdl4s.cwl

import eu.timepit.refined._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.string._
import shapeless.{:+:, CNil}

case class WorkflowStepInput(src: String)

case class InputParameter(
                           id: Option[String], //not really optional but can be specified upstream
                           label: Option[String],
                           secondaryFiles: Option[ECMAScriptExpression :+: String :+: Array[ECMAScriptExpression :+: String :+: CNil] :+: CNil],
                           format: Option[ECMAScriptExpression :+: String :+: Array[String] :+: CNil],
                           streamable: Option[Boolean],
                           doc: Option[String :+: Array[String] :+: CNil],
                           inputBinding: Option[CommandLineBinding],
                           default: Option[String], //can be of type "Any" which... sucks.
                           `type`: Option[MyriadInputType]) {

  type `type` = MyriadInputType
  type Id = String
}

case class InputRecordSchema(
  `type`: String Refined MatchesRegex[W.`"record"`.T],
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
  `type`: String Refined MatchesRegex[W.`"enum"`.T],
  label: Option[String],
  inputBinding: Option[CommandLineBinding])

case class InputArraySchema(
  items: MyriadInputType,
  `type`: String Refined MatchesRegex[W.`"array"`.T],
  label: Option[String],
  inputBinding: Option[CommandLineBinding])

case class CommandLineBinding(
                               loadContents: Option[Boolean],
                               position: Option[Int],
                               prefix: Option[String],
                               separate: Option[String],
                               itemSeparator: Option[String],
                               valueFrom: Option[ECMAScriptExpression :+: String :+: CNil], // could be "Expression" to be evaluated
                               shellQuote: Option[Boolean])

case class WorkflowOutputParameter(
                                    id: Option[String], //Really not optional but can be declared upstream
                                    label: Option[String],
                                    secondaryFiles: Option[ECMAScriptExpression :+: String :+: Array[ECMAScriptExpression :+: String :+: CNil] :+: CNil],
                                    format: Option[ECMAScriptExpression :+: String :+: Array[String] :+: CNil],
                                    streamable: Option[Boolean],
                                    doc: Option[String :+: Array[String] :+: CNil],
                                    `type`: Option[MyriadOutputType]) {

  type `type` = MyriadOutputType
  type Id = String
}

case class InputBinding(position: Int, prefix: String)

case class OutputRecordSchema(
  `type`: String Refined MatchesRegex[W.`"record"`.T],
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
  `type`: String Refined MatchesRegex[W.`"enum"`.T],
  label: Option[String],
  outputBinding: Option[CommandOutputBinding])

/** @see <a href="http://www.commonwl.org/v1.0/Workflow.html#CommandOutputBinding">CommandOutputBinding</a> */
case class CommandOutputBinding(
                                 glob: Option[ECMAScriptExpression :+: String :+: Array[String] :+: CNil],
                                 loadContents: Option[Boolean],
                                 outputEval: Option[ECMAScriptExpression :+: String :+: CNil])

case class OutputArraySchema(
  items: MyriadOutputType,
  `type`: String Refined MatchesRegex[W.`"array"`.T],
  label: Option[String],
  outputBinding: Option[CommandOutputBinding])


case class InlineJavascriptRequirement(
  `class`: String Refined MatchesRegex[W.`"InlineJavascriptRequirement"`.T],
  expressionLib: Option[Array[String]])

case class SchemaDefRequirement(
  `class`: String Refined MatchesRegex[W.`"SchemaDefRequirement"`.T],
  types: Array[InputRecordSchema :+: InputEnumSchema :+: InputArraySchema :+: CNil])

//There is a large potential for regex refinements on these string types
case class DockerRequirement(
  `class`: String Refined MatchesRegex[W.`"DockerRequirement"`.T],
  dockerPull: Option[String], //TODO Refine to match a docker image regex?
  dockerLoad: Option[String],
  dockerFile: Option[String],
  dockerImport: Option[String],
  dockerImageId: Option[String],
  dockerOutputDirectory: Option[String]
  )

case class SoftwareRequirement(
  `class`: String Refined MatchesRegex[W.`"SoftwareRequirement"`.T],
  packages: Array[SoftwarePackage] :+: Map[SoftwarePackage#Package, SoftwarePackage#Specs] :+: Map[SoftwarePackage#Package, SoftwarePackage] :+: CNil
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
  `class`: String Refined MatchesRegex[W.`"InitialWorkDirRequirement"`.T],
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
    CNil

  )

/** 
 *  Short for "Directory Entry"
 *  @see <a href="http://www.commonwl.org/v1.0/CommandLineTool.html#Dirent">Dirent Specification</a>
 */
case class Dirent(
                   entry: ECMAScriptExpression :+: String :+: CNil,
                   entryName: Option[ECMAScriptExpression :+: String :+: CNil],
                   writable: Option[Boolean]
  )

//TODO
//Figure out how to declare Any type
case class EnvVarRequirement(
  `class`: String Refined MatchesRegex[W.`"EnvVarRequirement"`.T],
  envDef: 
    Array[EnvironmentDef] :+:
    Map[EnvironmentDef#EnvName, EnvironmentDef#EnvValue] :+:
    Map[EnvironmentDef#EnvName, EnvironmentDef] :+:
    CNil)

case class EnvironmentDef(envName: String, envValue: ECMAScriptExpression :+: String :+: CNil) {
  type EnvName = String
  type EnvValue = String
}

case class ShellCommandRequirement(`class`: String Refined MatchesRegex[W.`"ShellCommandRequirement"`.T])

case class ResourceRequirement(
                                `class`: String Refined MatchesRegex[W.`"ResourceRequirement"`.T],
                                coresMin: Long :+: ECMAScriptExpression :+: String :+: CNil,
                                coresMax: Int :+: ECMAScriptExpression :+: String :+: CNil,
                                ramMin: Long :+: ECMAScriptExpression :+: String :+: CNil,
                                ramMax: Long :+: ECMAScriptExpression :+: String :+: CNil,
                                tmpdirMin: Long :+: ECMAScriptExpression :+: String :+: CNil,
                                tmpdirMax: Long :+: ECMAScriptExpression :+: String :+: CNil,
                                outdirMin: Long :+: ECMAScriptExpression :+: String :+: CNil,
                                outdirMax: Long :+: ECMAScriptExpression :+: String :+: CNil
  )

case class SubworkflowFeatureRequirement(
  `class`: String Refined MatchesRegex[W.`"SubworkflowFeatureRequirement"`.T]
  )

case class ScatterFeatureRequirement(
  `class`: String Refined MatchesRegex[W.`"ScatterFeatureRequirement"`.T]
  )

case class MultipleInputFeatureRequirement(
  `class`: String Refined MatchesRegex[W.`"MultipleInputFeatureRequirement"`.T]
  )

case class StepInputExpressionRequirement(
  `class`: String Refined MatchesRegex[W.`"StepInputExpressionRequirement"`.T]
  )
