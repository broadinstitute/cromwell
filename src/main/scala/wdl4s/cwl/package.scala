package wdl4s

import eu.timepit.refined.api.Refined
import eu.timepit.refined.string._
import eu.timepit.refined._
import shapeless.{:+:, CNil}
import io.circe.Decoder

package object cwl {

  import CwlType._
  import CwlVersion._
  import ScatterMethod._

  implicit val cwlTypeDecoder = Decoder.enumDecoder(CwlType)
  implicit val cwlVersionDecoder = Decoder.enumDecoder(CwlVersion)
  implicit val scatterMethodDecoder = Decoder.enumDecoder(ScatterMethod)

  type WorkflowStepInputId = String

  type WorkflowStepInputSource = String :+: Array[String] :+: CNil

  //These are supposed to be valid ECMAScript Expressions.  See http://www.commonwl.org/v1.0/Workflow.html#Expressions
  type ECMAScriptExpression = String Refined MatchesRegex[W.`"$({.*}|{.*})"`.T]

  type Requirement = 
    InlineJavascriptRequirement :+:
    SchemaDefRequirement :+:
    DockerRequirement :+:
    SoftwareRequirement :+:
    InitialWorkDirRequirement :+:
    EnvVarRequirement :+:
    ShellCommandRequirement :+:
    ResourceRequirement :+:
    SubworkflowFeatureRequirement :+:
    ScatterFeatureRequirement :+:
    MultipleInputFeatureRequirement :+:
    StepInputExpressionRequirement :+:
    CNil

  type MyriadInputType = 
    CwlType :+:
    InputRecordSchema :+:
    InputEnumSchema :+:
    InputArraySchema :+:
    String :+:
    Array[
      CwlType :+:
      InputRecordSchema :+:
      InputEnumSchema :+:
      InputArraySchema :+:
      String :+:
      CNil
    ] :+:
    CNil

  type MyriadOutputType = 
    CwlType :+:
    OutputRecordSchema :+:
    OutputEnumSchema :+:
    OutputArraySchema :+:
    String :+:
    Array[
      CwlType :+:
      OutputRecordSchema :+:
      OutputEnumSchema :+:
      OutputArraySchema :+:
      String :+:
      CNil
    ] :+:
    CNil

  type MyriadCommandInputType = 
    CwlType :+:
    CommandInputRecordSchema :+:
    CommandInputEnumSchema :+:
    CommandInputArraySchema :+:
    String :+:
    Array[
      CwlType  :+:
      CommandInputRecordSchema :+:
      CommandInputEnumSchema :+:
      CommandInputArraySchema :+:
      String :+:
      CNil 
      ] :+:
    CNil 

}
