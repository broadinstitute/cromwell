package cwl

import cwl.CommandLineTool._
import shapeless.{:+:, CNil}
import cwl.CwlType.CwlType
import io.circe.Json
import spray.json.JsValue
import wom.types.WomType
import spray.json._

trait TypeAliases {

  type CwlAny =
    File :+:
      Json :+:
      CNil

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

  type Hint =
    Requirement :+:
      CwlAny :+:
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
      MyriadInputInnerType :+:
      Array[MyriadInputInnerType] :+:
      CNil

  type MyriadInputInnerType =
    CwlType :+:
      InputRecordSchema :+:
      InputEnumSchema :+:
      InputArraySchema :+:
      String :+:
      CNil

  type MyriadOutputType =
      MyriadOutputInnerType :+:
      Array[MyriadOutputInnerType] :+:
      CNil

  type MyriadOutputInnerType =
    CwlType :+:
      OutputRecordSchema :+:
      OutputEnumSchema :+:
      OutputArraySchema :+:
      String :+:
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
  
  type ScatterVariables = Option[String :+: Array[String] :+: CNil]
  
  type ResourceRequirementType = Long :+: Expression :+: String :+: CNil
}

object CwlAny {

  implicit class EnhancedJson(val json: Json) extends AnyVal {
    // Yes, urgh! This is disgusting but until we rewrite the entire JS coercion set for circe...
    def sprayJsonRepresentation: JsValue = json.noSpaces.parseJson
  }

  object Json { def unapply(cwlAny: CwlAny): Option[Json] = cwlAny.select[Json] }
  object File { def unapply(cwlAny: CwlAny): Option[File] = cwlAny.select[File] }
}

object MyriadInputType {
  object CwlType {
    def unapply(m: MyriadInputType): Option[CwlType] = {
      m.select[MyriadInputInnerType].flatMap(_.select[CwlType])
    }
  }

  object WomType {
    def unapply(m: MyriadInputType): Option[WomType] = m match {
      case CwlType(c) => Option(cwl.cwlTypeToWomType(c))
      case _ => None
    }
  }
}
