package wdl4s

import wdl4s.cwl._

import io.circe.syntax._
import io.circe._
import io.circe.parser._
import io.circe.shapes._
import io.circe.generic.auto._
import io.circe.generic.semiauto._
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
 * decodeCwl(firstTool) //returns Either[Error, Cwl]
 * }}}
 *
 *
 * It makes heavy use of Circe YAML/Json auto derivation feature and
 * Circe modules that support the Scala libraries shapeless and Refined.
 *
 * The [[https://oss.sonatype.org/service/local/repositories/releases/archive/com/chuusai/shapeless_2.12/2.3.2/shapeless_2.12-2.3.2-javadoc.jar/!/shapeless/Coproduct.html shapeless.coproduct]] feature allows us to specify a large
 * number of potential types to be parsed.  A.k.a. a "disjunction" or "OR" relationship amongst these types.
 *
 * The [[https://github.com/fthomas/refined/blob/master/modules/core/shared/src/main/scala/eu/timepit/refined/string.scala MatchesRegex]] "refined type" is used
 * to enforce structure upon String values in the CWL Yaml.  Shapeless' Witness type
 * is used to declare a type containing a String literal.
 *
 * @see <a href="http://www.commonwl.org/v1.0/">CWL Specification</a>
 * @see <a href="https://github.com/circe/circe">circe</a>
 * @see <a href="https://github.com/circe/circe-yaml">circe-yaml</a>
 * @see <a href="https://github.com/fthomas/refined">Refined</a>
 * @see <a href="https://github.com/milessabin/shapeless">Shapeless</a>
 */
package object cwl {

  type CirceError[A] = Either[io.circe.Error, A]
  type CirceRead[A] = Kleisli[CirceError, String, A]

  import CwlType._
  import CwlVersion._
  import ScatterMethod._

  type Yaml = String

  implicit val cwlTypeDecoder = Decoder.enumDecoder(CwlType)
  implicit val cwlVersionDecoder = Decoder.enumDecoder(CwlVersion)
  implicit val scatterMethodDecoder = Decoder.enumDecoder(ScatterMethod)
  /*
  implicit def encodeAdtNoDiscr[A, Repr <: Coproduct](implicit
    gen: Generic.Aux[A, Repr],
    encodeRepr: Encoder[Repr]
  ): Encoder[A] = encodeRepr.contramap(gen.to)

  implicit def decodeAdtNoDiscr[A, Repr <: Coproduct](implicit
    gen: Generic.Aux[A, Repr],
    decodeRepr: Decoder[Repr]
): Decoder[A] = decodeRepr.map(gen.from)
        */

  def decodeCwl: Yaml => Either[Error, Cwl] =
      YamlParser.
        parse(_).
        map(_.noSpaces).
        flatMap(json => decode[CommandLineTool](json) orElse decode[Workflow](json))

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

    type WorkflowInput =
      Map[InputParameter#Id, InputParameter] :+:
      Map[InputParameter#Id, InputParameter#`type`] :+:
      Array[InputParameter] :+:
      CNil

    type WorkflowOutput =
      Map[WorkflowOutputParameter#Id, WorkflowOutputParameter] :+:
      Map[WorkflowOutputParameter#Id, WorkflowOutputParameter#`type`] :+:
      Array[WorkflowOutputParameter] :+:
      CNil

    type WorkflowSteps =
      Map[String, WorkflowStep] :+:
      Array[WorkflowStep] :+:
      CNil
}
