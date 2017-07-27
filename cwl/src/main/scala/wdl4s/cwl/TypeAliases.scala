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
}
