package wdl4s

import wdl4s.cwl._

import io.circe.syntax._
import io.circe._
import io.circe.parser._
import io.circe.shapes._
import io.circe.generic.auto._
import io.circe.yaml.{parser => YamlParser}
import io.circe.Json

import io.circe.syntax._
import io.circe._
import io.circe.parser._
import io.circe.shapes._
import io.circe.generic.auto._
import shapeless._, poly._//, ops.union._, union._
import shapeless.ops.coproduct._
import cats._, implicits._//, instances._
import cats.data.Kleisli
import io.circe._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.string._
import eu.timepit.refined._

import io.circe.refined._
import io.circe.generic.semiauto._
/**
 * This package is intended to parse all CWL files.
 *
 * =Usage=
 * {{{
 * import wdl4s.cwl._
 *
 * val firstTool = """
 * cwlVersion: v1.0
 * class: CommandLineTool
 * baseCommand: echo
 * inputs:
 *   message:
 *     type: string
 *     inputBinding:
 *       position: 1
 * outputs: []
 * """
 * decodeCommandLine(firstTool) //returns Either[Error, CommandLineTool]
 * }}}
 *
 *
 * It makes heavy use of Circe YAML/Json auto derivation feature and
 * Circe modules that support the Scala libraries shapeless and Refined.
 *
 * The [[https://oss.sonatype.org/service/local/repositories/releases/archive/com/chuusai/shapeless_2.12/2.3.2/shapeless_2.12-2.3.2-javadoc.jar/!/shapeless/Coproduct.html shapeless.coproduct]] feature allows us to specify a large
 * number of potential types against which a
 *
 * @see <a href="">CWL Specification</a>
 * @see <a href="">circe</a>
 * @see <a href="">circe-yaml</a>
 * @see <a href="">Refined</a>
 * @see <a href="">Shapeless</a>
 */
package object cwl {

  type CirceError[A] = Either[io.circe.Error, A]
  type CirceRead[A] = Kleisli[CirceError, String, A]

  import CwlType._
  import CwlVersion._
  import ScatterMethod._


  implicit val cwlTypeDecoder = Decoder.enumDecoder(CwlType)
  implicit val cwlVersionDecoder = Decoder.enumDecoder(CwlVersion)
  implicit val scatterMethodDecoder = Decoder.enumDecoder(ScatterMethod)

  implicit private val workflowDecoder: Decoder[Workflow] = deriveDecoder[Workflow]

  implicit private val commandLineDecoder: Decoder[CommandLineTool] = deriveDecoder[CommandLineTool]

  private def yamlToJson: CirceRead[String] = Kleisli[CirceError, String, String](YamlParser.parse(_).map(_.noSpaces))

  private def jsonToWorkflow: CirceRead[Workflow] =
    Kleisli[CirceError, String, Workflow](decode[Workflow])

  private def jsonToCommandLineTool: CirceRead[CommandLineTool] =
    Kleisli[CirceError, String, CommandLineTool](decode[CommandLineTool])

  //All the Kleisli noise buys us these nice derived functions:
  def commandLineFromYamlKleisli: CirceRead[CommandLineTool] = yamlToJson andThen jsonToCommandLineTool

  def commandLineFromYaml: String => Either[Error, CommandLineTool] = commandLineFromYamlKleisli.run

  def workflowFromYamlKleisli: CirceRead[Workflow] = yamlToJson andThen jsonToWorkflow

  def workflowFromYaml: String => Either[Error, Workflow] = workflowFromYamlKleisli.run


  type WorkflowStepInputId = String

  type WorkflowStepInputSource = String :+: Array[String] :+: CNil

  /**
   * These are supposed to be valid ECMAScript Expressions.
   * See http://www.commonwl.org/v1.0/Workflow.html#Expressions
   */
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
