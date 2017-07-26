package wdl4s.cwl

import shapeless.{:+:, CNil}
import wdl4s.cwl.CwlType.CwlType
import eu.timepit.refined._

trait TypeAliases {

  type CwlAny = String

  type WorkflowStepInputId = String

  type WorkflowStepInputSource = String :+: Array[String] :+: CNil

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

  type WorkflowInput =
    Map[InputParameter#Id, InputParameter] :+:
      Map[InputParameter#Id, InputParameter#`type`] :+:
      Array[InputParameter] :+:
      CNil

  type WorkflowOutput =
    Map[WorkflowOutputParameter#Id, WorkflowOutputParameter] :+:
      Map[WorkflowOutputParameter#Id, MyriadOutputType] :+:
      Array[WorkflowOutputParameter] :+:
      CNil

  type WorkflowSteps =
    Map[String, WorkflowStep] :+:
      Array[WorkflowStep] :+:
      CNil

  type EVR = EnvVarRequirement.`class`.type  => EnvVarRequirement
  type IJR = W.`"InlineJavascriptRequirement"`.T => InlineJavascriptRequirement
  type SR = W.`"SoftwareRequirement"`.T => SoftwareRequirement
  type SWFR = W.`"SubworkflowFeatureRequirement"`.T => SubworkflowFeatureRequirement
  type SDR = W.`"SchemaDefRequirement"`.T => SchemaDefRequirement
  type DR = W.`"DockerRequirement"`.T => DockerRequirement
  type IWDR = W.`"InitialWorkDirRequirement"`.T => InitialWorkDirRequirement
  type SCR = W.`"ShellCommandRequirement"`.T => ShellCommandRequirement
  type RR = W.`"ResourceRequirement"`.T => ResourceRequirement
  type SFR = W.`"ScatterFeatureRequirement"`.T => ScatterFeatureRequirement
  type MIFR = W.`"MultipleInputFeatureRequirement"`.T => MultipleInputFeatureRequirement
  type SIER = W.`"StepInputExpressionRequirement"`.T => StepInputExpressionRequirement

  type Target =
    EVR :+:
    IJR :+:
    SR :+:
    SWFR :+:
    SDR :+:
    DR :+:
    IWDR :+:
    SCR :+:
    RR :+:
    SFR :+:
    MIFR :+:
    SIER :+:
    CNil

}
