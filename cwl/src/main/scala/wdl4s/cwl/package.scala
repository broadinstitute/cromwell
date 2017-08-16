package wdl4s

import eu.timepit.refined.api.Refined
import eu.timepit.refined.string._
import eu.timepit.refined._
import wdl4s.cwl.CwlType.CwlType
import wdl4s.cwl.CwlType._
import wdl4s.wdl.types._

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
package object cwl extends TypeAliases {

  def cwlTypeToWdlType : CwlType => WdlType = {
    case Null => WdlNothingType
    case Boolean => WdlBooleanType
    case Int => WdlIntegerType
    case Long => WdlIntegerType
    case Float => WdlFloatType
    case Double => WdlFloatType
    case String => WdlStringType
    case CwlType.File => WdlFileType
    case CwlType.Directory => ???
  }


  /**
    *
    * These are supposed to be valid ECMAScript Expressions.
    * See http://www.commonwl.org/v1.0/Workflow.html#Expressions
    */
  type ECMAScriptExpression = String Refined MatchesRegex[W.`"$([^)]*)"`.T]

  type WdlTypeMap = Map[String, WdlType]

}
