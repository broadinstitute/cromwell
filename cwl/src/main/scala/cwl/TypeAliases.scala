package cwl

import cwl.CwlType.CwlType
import cwl.ExpressionEvaluator.{ECMAScriptExpression, ECMAScriptFunction}
import io.circe.Json
import shapeless.{:+:, CNil}
import wom.types.WomType

trait TypeAliases {

  type Expression = ECMAScriptExpression :+: ECMAScriptFunction :+: CNil

  type StringOrExpression = Expression :+: String :+: CNil

  type CwlAny =
    FileOrDirectory :+:
      Array[FileOrDirectory] :+:
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

  type ResourceRequirementType = Long :+: Expression :+: String :+: CNil

  type SingleOrArrayOfStrings = String :+: Array[String] :+: CNil

  type ScatterVariables = Option[SingleOrArrayOfStrings]

  type FileOrDirectory = File :+: Directory :+: CNil

  type Glob = StringOrExpression :+: Array[String] :+: CNil

  type SecondaryFiles = StringOrExpression :+: Array[StringOrExpression] :+: CNil
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
