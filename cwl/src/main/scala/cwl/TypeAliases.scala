package cwl

import cwl.CommandLineTool._
import shapeless.{:+:, CNil}
import cwl.CwlType.CwlType
import wom.types.WomType

trait TypeAliases {

  type CwlAny = String

  type WorkflowStepInputId = String

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

  // TODO WOM: Record Schema as well as Directories are not included because they're not supported yet, although they should eventually be.
  // Removing them saves some compile time when building decoders for this type (in CwlInputParsing)
  type MyriadInputValuePrimitives = 
    String :+:
    Int :+:
    Long :+:
    File :+:
    Float :+:
    Double :+:
    Boolean :+:
    CNil

  type MyriadInputValue =
    MyriadInputValuePrimitives :+:
      Array[MyriadInputValuePrimitives] :+:
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

  object MyriadInputType {
    object CwlType {
      def unapply(m: MyriadInputType): Option[CwlType] = {
        m.select[CwlType]
      }
    }

    object WomType {
      def unapply(m: MyriadInputType): Option[WomType] = m match {
        case CwlType(c) => Some(cwl.cwlTypeToWdlType(c))
        case _ => None
      }
    }
  }

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
