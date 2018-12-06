package cwl

import cats.data.NonEmptyList
import cats.syntax.either._
import cats.syntax.show._
import common.Checked
import common.validation.Checked._
import cwl.CommandLineTool.{CommandInputParameter, CommandOutputParameter}
import cwl.CwlType.CwlType
import cwl.CwlVersion.CwlVersion
import cwl.ExpressionTool.{ExpressionToolInputParameter, ExpressionToolOutputParameter}
import cwl.InitialWorkDirRequirement.IwdrListingArrayEntry
import cwl.LinkMergeMethod.LinkMergeMethod
import cwl.ScatterMethod.ScatterMethod
import cwl.SchemaDefRequirement.SchemaDefTypes
import cwl.Workflow.{WorkflowInputParameter, WorkflowOutputParameter}
import io.circe.DecodingFailure._
import io.circe.Json._
import io.circe._
import io.circe.generic.semiauto._
import io.circe.literal._
import io.circe.refined._
import io.circe.shapes._
import shapeless.Coproduct

object CwlCodecs {
  /*
  Semi-auto codecs for the Cwl types
  https://circe.github.io/circe/codecs/semiauto-derivation.html

  Some types, such as BaseCommand are not listed, as they are actually one-for-one type aliases. Having two exact
  implicits for the same type will confuse the compiler, leading to cryptic error messages elsewhere.

  For example, a duplicated decoder such as:

  ```
  implicit lazy val decodeSingleOrArrayOfStrings: Decoder[SingleOrArrayOfStrings] = decodeCCons
  implicit lazy val decodeBaseCommand: Decoder[BaseCommand] = decodeCCons
  ```

  causes cryptic errors such as:

  ```
  [error] cwl/src/main/scala/cwl/CwlCodecs.scala:123:456: could not find Lazy implicit value of type io.circe.generic.decoding.DerivedDecoder[A]
  [error]   implicit lazy val decodeCommandInputParameter: Decoder[CommandInputParameter] = deriveDecoder

  (etc. etc. etc.)
  ```
   */
  // Semi-Auto decoders
  implicit lazy val decodeArgumentCommandLineBinding: Decoder[ArgumentCommandLineBinding] = deriveDecoder
  implicit lazy val decodeCommandInputParameter: Decoder[CommandInputParameter] = deriveDecoder
  implicit lazy val decodeCommandLineTool: Decoder[CommandLineTool] = deriveDecoder
  implicit lazy val decodeCommandOutputBinding: Decoder[CommandOutputBinding] = deriveDecoder
  implicit lazy val decodeCommandOutputParameter: Decoder[CommandOutputParameter] = deriveDecoder
  implicit lazy val decodeCwlAny: Decoder[CwlAny] = decodeCCons
  implicit lazy val decodeCwlType: Decoder[CwlType] = Decoder.enumDecoder(CwlType)
  implicit lazy val decodeCwlVersion: Decoder[CwlVersion] = Decoder.enumDecoder(CwlVersion)
  implicit lazy val decodeDirectory: Decoder[Directory] = deriveDecoder
  implicit lazy val decodeDockerRequirement: Decoder[DockerRequirement] = deriveDecoder
  implicit lazy val decodeEnvVarRequirement: Decoder[EnvVarRequirement] = deriveDecoder
  implicit lazy val decodeEnvironmentDef: Decoder[EnvironmentDef] = deriveDecoder
  implicit lazy val decodeExpression: Decoder[Expression] = decodeCCons
  implicit lazy val decodeExpressionDirent: Decoder[ExpressionDirent] = deriveDecoder
  implicit lazy val decodeExpressionTool: Decoder[ExpressionTool] = deriveDecoder
  implicit lazy val decodeExpressionToolInputParameter: Decoder[ExpressionToolInputParameter] = deriveDecoder
  implicit lazy val decodeExpressionToolOutputParameter: Decoder[ExpressionToolOutputParameter] = deriveDecoder
  implicit lazy val decodeFile: Decoder[File] = deriveDecoder
  implicit lazy val decodeFileOrDirectory: Decoder[FileOrDirectory] = decodeCCons
  implicit lazy val decodeGlob: Decoder[Glob] = decodeCCons
  implicit lazy val decodeHint: Decoder[Hint] = decodeCCons
  implicit lazy val decodeInitialWorkDirRequirement: Decoder[InitialWorkDirRequirement] = deriveDecoder
  implicit lazy val decodeInlineJavascriptRequirement: Decoder[InlineJavascriptRequirement] = deriveDecoder
  implicit lazy val decodeInputArraySchema: Decoder[InputArraySchema] = deriveDecoder
  implicit lazy val decodeInputBinding: Decoder[InputBinding] = deriveDecoder
  implicit lazy val decodeInputCommandLineBinding: Decoder[InputCommandLineBinding] = deriveDecoder
  implicit lazy val decodeInputEnumSchema: Decoder[InputEnumSchema] = deriveDecoder
  implicit lazy val decodeInputRecordField: Decoder[InputRecordField] = deriveDecoder
  implicit lazy val decodeInputRecordSchema: Decoder[InputRecordSchema] = deriveDecoder
  implicit lazy val decodeInputResourceRequirement: Decoder[DnaNexusInputResourceRequirement] = deriveDecoder
  implicit lazy val decodeIwdrListingArrayEntry: Decoder[IwdrListingArrayEntry] = decodeCCons
  implicit lazy val decodeLinkMergeMethod: Decoder[LinkMergeMethod] = Decoder.enumDecoder(LinkMergeMethod)
  implicit lazy val decodeMultipleInputFeatureRequirement: Decoder[MultipleInputFeatureRequirement] = deriveDecoder
  implicit lazy val decodeMyriadInputInnerType: Decoder[MyriadInputInnerType] = decodeCCons
  implicit lazy val decodeMyriadInputType: Decoder[MyriadInputType] = decodeCCons
  implicit lazy val decodeMyriadOutputType: Decoder[MyriadOutputType] = decodeCCons
  implicit lazy val decodeOutputArraySchema: Decoder[OutputArraySchema] = deriveDecoder
  implicit lazy val decodeOutputEnumSchema: Decoder[OutputEnumSchema] = deriveDecoder
  implicit lazy val decodeOutputRecordField: Decoder[OutputRecordField] = deriveDecoder
  implicit lazy val decodeOutputRecordSchema: Decoder[OutputRecordSchema] = deriveDecoder
  implicit lazy val decodeRequirement: Decoder[Requirement] = decodeCCons
  implicit lazy val decodeResourceRequirement: Decoder[ResourceRequirement] = deriveDecoder
  implicit lazy val decodeResourceRequirementType: Decoder[ResourceRequirementType] = decodeCCons
  implicit lazy val decodeScatterFeatureRequirement: Decoder[ScatterFeatureRequirement] = deriveDecoder
  implicit lazy val decodeScatterMethod: Decoder[ScatterMethod] = Decoder.enumDecoder(ScatterMethod)
  implicit lazy val decodeSchemaDefRequirement: Decoder[SchemaDefRequirement] = deriveDecoder
  implicit lazy val decodeSchemaDefTypes: Decoder[SchemaDefTypes] = decodeCCons
  implicit lazy val decodeSecondaryFiles: Decoder[SecondaryFiles] = decodeCCons
  implicit lazy val decodeShellCommandRequirement: Decoder[ShellCommandRequirement] = deriveDecoder
  implicit lazy val decodeSingleOrArrayOfStrings: Decoder[SingleOrArrayOfStrings] = decodeCCons
  implicit lazy val decodeSoftwarePackage: Decoder[SoftwarePackage] = deriveDecoder
  implicit lazy val decodeSoftwareRequirement: Decoder[SoftwareRequirement] = deriveDecoder
  implicit lazy val decodeStepInputExpressionRequirement: Decoder[StepInputExpressionRequirement] = deriveDecoder
  implicit lazy val decodeStringDirent: Decoder[StringDirent] = deriveDecoder
  implicit lazy val decodeStringOrExpression: Decoder[StringOrExpression] = decodeCCons
  implicit lazy val decodeSubworkflowFeatureRequirement: Decoder[SubworkflowFeatureRequirement] = deriveDecoder
  implicit lazy val decodeWorkflow: Decoder[Workflow] = deriveDecoder
  implicit lazy val decodeWorkflowInputParameter: Decoder[WorkflowInputParameter] = deriveDecoder
  implicit lazy val decodeWorkflowOutputParameter: Decoder[WorkflowOutputParameter] = deriveDecoder
  implicit lazy val decodeWorkflowStep: Decoder[WorkflowStep] = deriveDecoder
  implicit lazy val decodeWorkflowStepInput: Decoder[WorkflowStepInput] = deriveDecoder
  implicit lazy val decodeWorkflowStepOutput: Decoder[WorkflowStepOutput] = deriveDecoder

  // Semi-Auto encoders
  implicit lazy val encodeArgumentCommandLineBinding: Encoder[ArgumentCommandLineBinding] = deriveEncoder
  implicit lazy val encodeCommandInputParameter: Encoder[CommandInputParameter] = deriveEncoder
  implicit lazy val encodeCommandLineTool: Encoder[CommandLineTool] = deriveEncoder
  implicit lazy val encodeCommandOutputBinding: Encoder[CommandOutputBinding] = deriveEncoder
  implicit lazy val encodeCommandOutputParameter: Encoder[CommandOutputParameter] = deriveEncoder
  implicit lazy val encodeCwlAny: Encoder[CwlAny] = encodeCCons
  implicit lazy val encodeCwlType: Encoder[CwlType] = Encoder.enumEncoder(CwlType)
  implicit lazy val encodeCwlVersion: Encoder[CwlVersion] = Encoder.enumEncoder(CwlVersion)
  implicit lazy val encodeDirectory: Encoder[Directory] = deriveEncoder
  implicit lazy val encodeDockerRequirement: Encoder[DockerRequirement] = deriveEncoder
  implicit lazy val encodeEnvVarRequirement: Encoder[EnvVarRequirement] = deriveEncoder
  implicit lazy val encodeEnvironmentDef: Encoder[EnvironmentDef] = deriveEncoder
  implicit lazy val encodeExpression: Encoder[Expression] = encodeCCons
  implicit lazy val encodeExpressionDirent: Encoder[ExpressionDirent] = deriveEncoder
  implicit lazy val encodeExpressionTool: Encoder[ExpressionTool] = deriveEncoder
  implicit lazy val encodeExpressionToolInputParameter: Encoder[ExpressionToolInputParameter] = deriveEncoder
  implicit lazy val encodeExpressionToolOutputParameter: Encoder[ExpressionToolOutputParameter] = deriveEncoder
  implicit lazy val encodeFile: Encoder[File] = deriveEncoder
  implicit lazy val encodeFileOrDirectory: Encoder[FileOrDirectory] = encodeCCons
  implicit lazy val encodeGlob: Encoder[Glob] = encodeCCons
  implicit lazy val encodeHint: Encoder[Hint] = encodeCCons
  implicit lazy val encodeInitialWorkDirRequirement: Encoder[InitialWorkDirRequirement] = deriveEncoder
  implicit lazy val encodeInlineJavascriptRequirement: Encoder[InlineJavascriptRequirement] = deriveEncoder
  implicit lazy val encodeInputArraySchema: Encoder[InputArraySchema] = deriveEncoder
  implicit lazy val encodeInputBinding: Encoder[InputBinding] = deriveEncoder
  implicit lazy val encodeInputCommandLineBinding: Encoder[InputCommandLineBinding] = deriveEncoder
  implicit lazy val encodeInputEnumSchema: Encoder[InputEnumSchema] = deriveEncoder
  implicit lazy val encodeInputRecordField: Encoder[InputRecordField] = deriveEncoder
  implicit lazy val encodeInputRecordSchema: Encoder[InputRecordSchema] = deriveEncoder
  implicit lazy val encodeInputResourceRequirement: Encoder[DnaNexusInputResourceRequirement] = deriveEncoder
  implicit lazy val encodeIwdrListingArrayEntry: Encoder[IwdrListingArrayEntry] = encodeCCons
  implicit lazy val encodeLinkMergeMethod: Encoder[LinkMergeMethod] = Encoder.enumEncoder(LinkMergeMethod)
  implicit lazy val encodeMultipleInputFeatureRequirement: Encoder[MultipleInputFeatureRequirement] = deriveEncoder
  implicit lazy val encodeMyriadInputInnerType: Encoder[MyriadInputInnerType] = encodeCCons
  implicit lazy val encodeMyriadInputType: Encoder[MyriadInputType] = encodeCCons
  implicit lazy val encodeMyriadOutputType: Encoder[MyriadOutputType] = encodeCCons
  implicit lazy val encodeOutputArraySchema: Encoder[OutputArraySchema] = deriveEncoder
  implicit lazy val encodeOutputEnumSchema: Encoder[OutputEnumSchema] = deriveEncoder
  implicit lazy val encodeOutputRecordField: Encoder[OutputRecordField] = deriveEncoder
  implicit lazy val encodeOutputRecordSchema: Encoder[OutputRecordSchema] = deriveEncoder
  implicit lazy val encodeRequirement: Encoder[Requirement] = encodeCCons
  implicit lazy val encodeResourceRequirement: Encoder[ResourceRequirement] = deriveEncoder
  implicit lazy val encodeResourceRequirementType: Encoder[ResourceRequirementType] = encodeCCons
  implicit lazy val encodeScatterFeatureRequirement: Encoder[ScatterFeatureRequirement] = deriveEncoder
  implicit lazy val encodeScatterMethod: Encoder[ScatterMethod] = Encoder.enumEncoder(ScatterMethod)
  implicit lazy val encodeSchemaDefRequirement: Encoder[SchemaDefRequirement] = deriveEncoder
  implicit lazy val encodeSchemaDefTypes: Encoder[SchemaDefTypes] = encodeCCons
  implicit lazy val encodeSecondaryFiles: Encoder[SecondaryFiles] = encodeCCons
  implicit lazy val encodeShellCommandRequirement: Encoder[ShellCommandRequirement] = deriveEncoder
  implicit lazy val encodeSingleOrArrayOfStrings: Encoder[SingleOrArrayOfStrings] = encodeCCons
  implicit lazy val encodeSoftwarePackage: Encoder[SoftwarePackage] = deriveEncoder
  implicit lazy val encodeSoftwareRequirement: Encoder[SoftwareRequirement] = deriveEncoder
  implicit lazy val encodeStepInputExpressionRequirement: Encoder[StepInputExpressionRequirement] = deriveEncoder
  implicit lazy val encodeStringDirent: Encoder[StringDirent] = deriveEncoder
  implicit lazy val encodeStringOrExpression: Encoder[StringOrExpression] = encodeCCons
  implicit lazy val encodeSubworkflowFeatureRequirement: Encoder[SubworkflowFeatureRequirement] = deriveEncoder
  implicit lazy val encodeWorkflow: Encoder[Workflow] = deriveEncoder
  implicit lazy val encodeWorkflowInputParameter: Encoder[WorkflowInputParameter] = deriveEncoder
  implicit lazy val encodeWorkflowOutputParameter: Encoder[WorkflowOutputParameter] = deriveEncoder
  implicit lazy val encodeWorkflowStep: Encoder[WorkflowStep] = deriveEncoder
  implicit lazy val encodeWorkflowStepInput: Encoder[WorkflowStepInput] = deriveEncoder
  implicit lazy val encodeWorkflowStepOutput: Encoder[WorkflowStepOutput] = deriveEncoder

  def decodeCwl(json: Json): Checked[Cwl] = {
    findClass(json) match {
      case Some("Workflow") => decodeWithErrorStringsJson[Workflow](json).map(Coproduct[Cwl].apply(_))
      case Some("CommandLineTool") => decodeWithErrorStringsJson[CommandLineTool](json).map(Coproduct[Cwl].apply(_))
      case Some("ExpressionTool") => decodeWithErrorStringsJson[ExpressionTool](json).map(Coproduct[Cwl].apply(_))
      case Some(other) => s"Class field was declared incorrectly: $other is not one of Workflow, CommandLineTool, or ExpressionTool! as seen in ${json.show}".invalidNelCheck
      case None => s"Class field was omitted in ${json.show}".invalidNelCheck
    }
  }

  private def decodeWithErrorStringsJson[A](in: Json)(implicit d: Decoder[A]): Checked[A] =
    in.as[A].leftMap(_.show).leftMap(NonEmptyList.one)

  private def findClass(json: Json): Option[String] =
    for {
      obj <- json.asObject
      map = obj.toMap
      classObj <- map.get("class")
      classString <- classObj.asString
    } yield classString
}
