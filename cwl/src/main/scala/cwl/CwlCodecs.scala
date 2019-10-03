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
  implicit lazy val decodeArgumentCommandLineBinding: Codec[ArgumentCommandLineBinding] = deriveCodec
  implicit lazy val decodeCommandInputParameter: Codec[CommandInputParameter] = deriveCodec
  implicit lazy val decodeCommandLineTool: Codec[CommandLineTool] = deriveCodec
  implicit lazy val decodeCommandOutputBinding: Codec[CommandOutputBinding] = deriveCodec
  implicit lazy val decodeCommandOutputParameter: Codec[CommandOutputParameter] = deriveCodec
  implicit lazy val decodeCwlType: Codec[CwlType] = Codec.codecForEnumeration(CwlType)
  implicit lazy val decodeCwlVersion: Codec[CwlVersion] = Codec.codecForEnumeration(CwlVersion)
  implicit lazy val decodeDirectory: Codec[Directory] = deriveCodec
  implicit lazy val decodeDockerRequirement: Codec[DockerRequirement] = deriveCodec
  implicit lazy val decodeEnvVarRequirement: Codec[EnvVarRequirement] = deriveCodec
  implicit lazy val decodeEnvironmentDef: Codec[EnvironmentDef] = deriveCodec
  implicit lazy val decodeExpressionDirent: Codec[ExpressionDirent] = deriveCodec
  implicit lazy val decodeExpressionTool: Codec[ExpressionTool] = deriveCodec
  implicit lazy val decodeExpressionToolInputParameter: Codec[ExpressionToolInputParameter] = deriveCodec
  implicit lazy val decodeExpressionToolOutputParameter: Codec[ExpressionToolOutputParameter] = deriveCodec
  implicit lazy val decodeFile: Codec[File] = deriveCodec
  implicit lazy val decodeInitialWorkDirRequirement: Codec[InitialWorkDirRequirement] = deriveCodec
  implicit lazy val decodeInlineJavascriptRequirement: Codec[InlineJavascriptRequirement] = deriveCodec
  implicit lazy val decodeInputArraySchema: Codec[InputArraySchema] = deriveCodec
  implicit lazy val decodeInputBinding: Codec[InputBinding] = deriveCodec
  implicit lazy val decodeInputCommandLineBinding: Codec[InputCommandLineBinding] = deriveCodec
  implicit lazy val decodeInputEnumSchema: Codec[InputEnumSchema] = deriveCodec
  implicit lazy val decodeInputRecordField: Codec[InputRecordField] = deriveCodec
  implicit lazy val decodeInputRecordSchema: Codec[InputRecordSchema] = deriveCodec
  implicit lazy val decodeInputResourceRequirement: Codec[DnaNexusInputResourceRequirement] = deriveCodec
  implicit lazy val decodeLinkMergeMethod: Codec[LinkMergeMethod] = Codec.codecForEnumeration(LinkMergeMethod)
  implicit lazy val decodeMultipleInputFeatureRequirement: Codec[MultipleInputFeatureRequirement] = deriveCodec
  implicit lazy val decodeOutputArraySchema: Codec[OutputArraySchema] = deriveCodec
  implicit lazy val decodeOutputEnumSchema: Codec[OutputEnumSchema] = deriveCodec
  implicit lazy val decodeOutputRecordField: Codec[OutputRecordField] = deriveCodec
  implicit lazy val decodeOutputRecordSchema: Codec[OutputRecordSchema] = deriveCodec
  implicit lazy val decodeResourceRequirement: Codec[ResourceRequirement] = deriveCodec
  implicit lazy val decodeScatterFeatureRequirement: Codec[ScatterFeatureRequirement] = deriveCodec
  implicit lazy val decodeScatterMethod: Codec[ScatterMethod] = Codec.codecForEnumeration(ScatterMethod)
  implicit lazy val decodeSchemaDefRequirement: Codec[SchemaDefRequirement] = deriveCodec
  implicit lazy val decodeShellCommandRequirement: Codec[ShellCommandRequirement] = deriveCodec
  implicit lazy val decodeSoftwarePackage: Codec[SoftwarePackage] = deriveCodec
  implicit lazy val decodeSoftwareRequirement: Codec[SoftwareRequirement] = deriveCodec
  implicit lazy val decodeStepInputExpressionRequirement: Codec[StepInputExpressionRequirement] = deriveCodec
  implicit lazy val decodeStringDirent: Codec[StringDirent] = deriveCodec
  implicit lazy val decodeSubworkflowFeatureRequirement: Codec[SubworkflowFeatureRequirement] = deriveCodec
  implicit lazy val decodeWorkflow: Codec[Workflow] = deriveCodec
  implicit lazy val decodeWorkflowInputParameter: Codec[WorkflowInputParameter] = deriveCodec
  implicit lazy val decodeWorkflowOutputParameter: Codec[WorkflowOutputParameter] = deriveCodec
  implicit lazy val decodeWorkflowStep: Codec[WorkflowStep] = deriveCodec
  implicit lazy val decodeWorkflowStepInput: Codec[WorkflowStepInput] = deriveCodec
  implicit lazy val decodeWorkflowStepOutput: Codec[WorkflowStepOutput] = deriveCodec

  // Semi-Auto encoders
  implicit lazy val decodeCwlAny: Decoder[CwlAny] = decodeCCons
  implicit lazy val decodeExpression: Decoder[Expression] = decodeCCons
  implicit lazy val decodeFileOrDirectory: Decoder[FileOrDirectory] = decodeCCons
  implicit lazy val decodeGlob: Decoder[Glob] = decodeCCons
  implicit lazy val decodeHint: Decoder[Hint] = decodeCCons
  implicit lazy val decodeIwdrListingArrayEntry: Decoder[IwdrListingArrayEntry] = decodeCCons
  implicit lazy val decodeMyriadInputInnerType: Decoder[MyriadInputInnerType] = decodeCCons
  implicit lazy val decodeMyriadInputType: Decoder[MyriadInputType] = decodeCCons
  implicit lazy val decodeMyriadOutputType: Decoder[MyriadOutputType] = decodeCCons
  implicit lazy val decodeRequirement: Decoder[Requirement] = decodeCCons
  implicit lazy val decodeResourceRequirementType: Decoder[ResourceRequirementType] = decodeCCons
  implicit lazy val decodeSchemaDefTypes: Decoder[SchemaDefTypes] = decodeCCons
  implicit lazy val decodeSecondaryFiles: Decoder[SecondaryFiles] = decodeCCons
  implicit lazy val decodeSingleOrArrayOfStrings: Decoder[SingleOrArrayOfStrings] = decodeCCons
  implicit lazy val decodeStringOrExpression: Decoder[StringOrExpression] = decodeCCons

  implicit lazy val encodeCwlAny: Encoder[CwlAny] = encodeCCons
  implicit lazy val encodeExpression: Encoder[Expression] = encodeCCons
  implicit lazy val encodeFileOrDirectory: Encoder[FileOrDirectory] = encodeCCons
  implicit lazy val encodeGlob: Encoder[Glob] = encodeCCons
  implicit lazy val encodeHint: Encoder[Hint] = encodeCCons
  implicit lazy val encodeIwdrListingArrayEntry: Encoder[IwdrListingArrayEntry] = encodeCCons
  implicit lazy val encodeMyriadInputInnerType: Encoder[MyriadInputInnerType] = encodeCCons
  implicit lazy val encodeMyriadInputType: Encoder[MyriadInputType] = encodeCCons
  implicit lazy val encodeMyriadOutputType: Encoder[MyriadOutputType] = encodeCCons
  implicit lazy val encodeRequirement: Encoder[Requirement] = encodeCCons
  implicit lazy val encodeResourceRequirementType: Encoder[ResourceRequirementType] = encodeCCons
  implicit lazy val encodeSchemaDefTypes: Encoder[SchemaDefTypes] = encodeCCons
  implicit lazy val encodeSecondaryFiles: Encoder[SecondaryFiles] = encodeCCons
  implicit lazy val encodeSingleOrArrayOfStrings: Encoder[SingleOrArrayOfStrings] = encodeCCons
  implicit lazy val encodeStringOrExpression: Encoder[StringOrExpression] = encodeCCons

  def decodeCwl(json: Json): Checked[Cwl] = {
    findClass(json) match {
      case Some("Workflow") => decodeWithErrorStringsJson[Workflow](json).map(Coproduct[Cwl].apply(_))
      case Some("CommandLineTool") => decodeWithErrorStringsJson[CommandLineTool](json).map(Coproduct[Cwl].apply(_))
      case Some("ExpressionTool") => decodeWithErrorStringsJson[ExpressionTool](json).map(Coproduct[Cwl].apply(_))
      case Some(other) => s"Class field was declared incorrectly: $other is not one of Workflow, CommandLineTool, or ExpressionTool! as seen in ${json.show}".invalidNelCheck
      case None => s"Class field was omitted in ${json.show}".invalidNelCheck
    }
  }

  private def decodeWithErrorStringsJson[A](in: Json)(implicit d: Codec[A]): Checked[A] =
    in.as[A].leftMap(_.show).leftMap(NonEmptyList.one)

  private def findClass(json: Json): Option[String] =
    for {
      obj <- json.asObject
      map = obj.toMap
      classObj <- map.get("class")
      classString <- classObj.asString
    } yield classString
}
