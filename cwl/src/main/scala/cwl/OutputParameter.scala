package cwl

import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import shapeless.{:+:, CNil, Poly1}
import wom.types.WomFileType
import wom.values.{WomFile, WomString, WomValue}

trait OutputParameter {
  def id: String
  def label: Option[String]
  def secondaryFiles: Option[SecondaryFiles]
  def format: Option[OutputParameterFormat]
  def streamable: Option[Boolean]
  def doc: Option[String :+: Array[String] :+: CNil]
  def outputBinding: Option[CommandOutputBinding]
  def `type`: Option[MyriadOutputType]
}

object OutputParameter {
  object IdAndType {
    def unapply(arg: OutputParameter): Option[(String, MyriadOutputType)] = arg.`type`.map((arg.id, _))
  }

  import cats.instances.list._
  import cats.instances.option._
  def format(formatOption: Option[StringOrExpression],
             parameterContext: ParameterContext,
             expressionLib: ExpressionLib): ErrorOr[Option[String]] = {
    formatOption.traverse[ErrorOr, String] {
      format(_, parameterContext, expressionLib)
    }
  }

  def format(format: StringOrExpression, parameterContext: ParameterContext, expressionLib: ExpressionLib): ErrorOr[String] = {
    format.fold(OutputParameter.FormatPoly).apply(parameterContext, expressionLib)
  }

  type FormatFunction = (ParameterContext, ExpressionLib) => ErrorOr[String]

  object FormatPoly extends Poly1 {
    implicit def caseStringOrExpression: Case.Aux[StringOrExpression, FormatFunction] = {
      at {
        _.fold(this)
      }
    }

    implicit def caseExpression: Case.Aux[Expression, FormatFunction] = {
      at {
        expression =>
          (parameterContext, expressionLib) =>
            val result: ErrorOr[WomValue] = ExpressionEvaluator.eval(expression, parameterContext, expressionLib)
            result flatMap {
              case womString: WomString => womString.value.valid
              case other => s"Not a valid file format: $other".invalidNel
            }
      }
    }

    implicit def caseString: Case.Aux[String, FormatFunction] = at { string => (_,_) => string.valid }
  }

  /**
    * Returns the list of secondary files for the primary file.
    */
  def secondaryFiles(primaryWomFile: WomFile,
                     stringWomFileType: WomFileType,
                     secondaryFilesOption: Option[SecondaryFiles],
                     parameterContext: ParameterContext,
                     expressionLib: ExpressionLib): ErrorOr[List[WomFile]] = {
    secondaryFilesOption
      .map(secondaryFiles(primaryWomFile, stringWomFileType, _, parameterContext, expressionLib))
      .getOrElse(Nil.valid)
  }

  /**
    * Returns the list of secondary files for the primary file.
    */
  def secondaryFiles(primaryWomFile: WomFile,
                     stringWomFileType: WomFileType,
                     secondaryFiles: SecondaryFiles,
                     parameterContext: ParameterContext,
                     expressionLib: ExpressionLib): ErrorOr[List[WomFile]] = {
    secondaryFiles
      .fold(OutputParameter.SecondaryFilesPoly)
      .apply(primaryWomFile, stringWomFileType, parameterContext, expressionLib)
  }

  type SecondaryFilesFunction = (WomFile, WomFileType, ParameterContext, ExpressionLib) => ErrorOr[List[WomFile]]

  object SecondaryFilesPoly extends Poly1 {
    implicit def caseStringOrExpression: Case.Aux[StringOrExpression, SecondaryFilesFunction] = {
      at {
        _.fold(this)
      }
    }

    implicit def caseExpression: Case.Aux[Expression, SecondaryFilesFunction] = {
      at {
        expression =>
          (primaryWomFile, stringWomFileType, parameterContext, expressionLib) =>
            File.secondaryExpressionFiles(primaryWomFile, stringWomFileType, expression, parameterContext, expressionLib)
      }
    }

    implicit def caseString: Case.Aux[String, SecondaryFilesFunction] = {
      at {
        string =>
          (primaryWomFile, stringWomFileType, _, _) =>
            File.secondaryStringFile(primaryWomFile, stringWomFileType, string).map(List(_))
      }
    }

    implicit def caseArray: Case.Aux[Array[StringOrExpression], SecondaryFilesFunction] = {
      at {
        array =>
          (primaryWomFile, stringWomFileType, parameterContext, expressionLib) =>
            val functions: List[SecondaryFilesFunction] = array.toList.map(_.fold(this))
            functions.flatTraverse(_ (primaryWomFile, stringWomFileType, parameterContext, expressionLib))
      }
    }
  }
}
