
import cwl.CwlType.{CwlType, _}
import cwl.ExpressionEvaluator.{ECMAScriptExpression, ECMAScriptFunction}
import lenthall.Checked
import shapeless._
import wom.executable.Executable
import wom.types._

/**
 * This package is intended to parse all CWL files.
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

  type Cwl = Workflow :+: CommandLineTool :+: CNil

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

  type StringOrExpression = Expression :+: String :+: CNil

  type Expression = ECMAScriptExpression :+: ECMAScriptFunction :+: CNil


  type WdlTypeMap = Map[String, WdlType]

  object CwlToWomExecutable extends Poly1 {
    implicit def caseClt = at[CommandLineTool](clt => clt.womExecutable())
    implicit def caseWf = at[Workflow](wf => wf.womExecutable())
  }

  implicit class CwlHelper(val cwl: Cwl) extends AnyVal {
    def womExecutable: Checked[Executable] = cwl.fold(CwlToWomExecutable)
  }
}
