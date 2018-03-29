
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr._
import cwl.CwlType._
import cwl.ExpressionEvaluator.{ECMAScriptExpression, ECMAScriptFunction}
import cwl.command.ParentName
import shapeless._
import wom.executable.Executable
import wom.expression.{IoFunctionSet, NoIoFunctionSet}
import wom.types._

import scala.util.{Failure, Success, Try}

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

  type CwlFile = Array[Cwl] :+: Cwl :+: CNil
  type Cwl = Workflow :+: CommandLineTool :+: ExpressionTool :+: CNil
  
  object Cwl {
    object Workflow { def unapply(cwl: Cwl): Option[Workflow] = cwl.select[Workflow] }
    object CommandLineTool { def unapply(cwl: Cwl): Option[CommandLineTool] = cwl.select[CommandLineTool] }
    object ExpressionTool { def unapply(cwl: Cwl): Option[ExpressionTool] = cwl.select[ExpressionTool] }
  }

  def cwlTypeToWomType : CwlType => WomType = {
    case CwlType.Any => WomAnyType
    case Null => WomNothingType
    case Boolean => WomBooleanType
    case Int => WomIntegerType
    case Long => WomLongType
    case Float => WomFloatType
    case Double => WomFloatType
    case String => WomStringType
    case CwlType.File => WomMaybePopulatedFileType
    case CwlType.Directory => WomMaybeListedDirectoryType
  }

  object StringOrExpression {
    object String {
      def unapply(soe: StringOrExpression): Option[String] = soe.select[String]
    }
    object Expression {
      def unapply(soe: StringOrExpression): Option[Expression] = soe.select[Expression]
    }
    object ECMAScriptExpression {
      def unapply(soe: StringOrExpression): Option[ECMAScriptExpression] = soe.select[Expression].flatMap(_.select[ECMAScriptExpression])
    }
    object ECMAScriptFunction {
      def unapply(soe: StringOrExpression): Option[ECMAScriptFunction] = soe.select[Expression].flatMap(_.select[ECMAScriptFunction])
    }
  }

  object Expression {
    object ECMAScriptExpression {
      def unapply(soe: Expression): Option[ECMAScriptExpression] = soe.select[ECMAScriptExpression]
    }
    object ECMAScriptFunction {
      def unapply(soe: Expression): Option[ECMAScriptFunction] = soe.select[ECMAScriptFunction]
    }
  }

  type WomTypeMap = Map[String, WomType]

  type RequirementsValidator = Requirement => ErrorOr[Requirement]
  import cats.syntax.validated._
  val AcceptAllRequirements: RequirementsValidator = _.validNel

  implicit class CwlHelper(val cwl: Cwl) extends AnyVal {
    def womExecutable(validator: RequirementsValidator, inputsFile: Option[String], ioFunctions: IoFunctionSet, strictValidation: Boolean): Checked[Executable] = {
      def executable = cwl match {
        case Cwl.Workflow(w) => w.womExecutable(validator, inputsFile, ioFunctions, strictValidation)
        case Cwl.CommandLineTool(clt) => clt.womExecutable(validator, inputsFile, ioFunctions, strictValidation)
        case Cwl.ExpressionTool(et) => et.womExecutable(validator, inputsFile, ioFunctions, strictValidation)
      }
      Try(executable) match {
        case Success(s) => s
        case Failure(f) => f.getMessage.invalidNelCheck
      }
    }

    def requiredInputs: Map[String, WomType] = {
      implicit val parent = ParentName.empty
      
      cwl match {
        case Cwl.Workflow(w) => selectWomTypeInputs(w.inputs collect {
          case i if i.`type`.isDefined => FullyQualifiedName(i.id).id -> i.`type`.get
        })
        case Cwl.CommandLineTool(clt) => selectWomTypeInputs(clt.inputs collect {
          case i if i.`type`.isDefined => FullyQualifiedName(i.id).id -> i.`type`.get
        })
      }
    }

    private def selectWomTypeInputs(myriadInputMap: Array[(String, MyriadInputType)]): Map[String, WomType] = {
      (myriadInputMap collect {
        case (key, MyriadInputType.WomType(w)) => key -> w
      }).toMap
    }
  }

  type ExpressionLib = Vector[String]
}
