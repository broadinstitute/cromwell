package cwl

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation.validate
import shapeless.Poly1
import wom.expression.IoFunctionSet
import wom.types.{WomFileType, WomSingleFileType}
import wom.values.{WomArray, WomFile, WomMaybePopulatedFile, WomValue}

import scala.concurrent.Await
import scala.concurrent.duration._

object FileParameter {
  private val ReadLimit = Option(64 * 1024)
  private val ReadTimeout = 60.seconds

  def populateSecondaryFiles(womValue: WomValue,
                             secondaryFilesCoproduct: Option[SecondaryFiles],
                             formatOption: Option[String],
                             parameterContext: ParameterContext,
                             expressionLib: ExpressionLib): ErrorOr[WomValue] = {

    womValue match {

      case womMaybePopulatedFile: WomMaybePopulatedFile =>
        val secondaryFilesErrorOr = FileParameter.secondaryFiles(
          womMaybePopulatedFile,
          WomSingleFileType,
          secondaryFilesCoproduct,
          parameterContext,
          expressionLib
        )

        secondaryFilesErrorOr map { secondaryFiles =>
          womMaybePopulatedFile.copy(secondaryFiles = secondaryFiles, formatOption = formatOption)
        }

      case womArray: WomArray =>
        womArray.value.toList.traverse(
          populateSecondaryFiles(_, secondaryFilesCoproduct, formatOption, parameterContext, expressionLib)
        ).map(WomArray(_))

      case womValue: WomValue => womValue.valid
    }
  }


  def load64KiB(path: String, ioFunctionSet: IoFunctionSet): ErrorOr[String] = {
    // TODO: WOM: propagate Future (or IO?) signature
    validate {
      val content = ioFunctionSet.readFile(path, ReadLimit, failOnOverflow = false)
      Await.result(content, ReadTimeout)
    }
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
      .fold(SecondaryFilesPoly)
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
