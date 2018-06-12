package cwl

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation.validate
import cwl.ontology.Schema
import shapeless.Poly1
import wom.expression.IoFunctionSet
import wom.types.{WomFileType, WomMaybePopulatedFileType}
import wom.values.{WomArray, WomFile, WomMaybePopulatedFile, WomValue}

import scala.concurrent.duration._
import scala.concurrent.{Await, Future}
import scala.util.Try

object FileParameter {
  private val ReadLimit = Option(64 * 1024)
  val ReadTimeout = 60.seconds

  def sync[A](f: Future[A]): Try[A] = Try(Await.result(f, FileParameter.ReadTimeout))

  def populateSecondaryFiles(womValue: WomValue,
                             secondaryFilesCoproduct: Option[SecondaryFiles],
                             formatOption: Option[String],
                             parameterContext: ParameterContext,
                             expressionLib: ExpressionLib,
                             ioFunctions: IoFunctionSet): ErrorOr[WomValue] = {

    womValue match {

      case womMaybePopulatedFile: WomMaybePopulatedFile =>
        val secondaryFilesErrorOr = FileParameter.secondaryFiles(
          womMaybePopulatedFile,
          WomMaybePopulatedFileType,
          secondaryFilesCoproduct,
          parameterContext,
          expressionLib,
          ioFunctions
        )

        secondaryFilesErrorOr map { secondaryFiles =>
          womMaybePopulatedFile.copy(secondaryFiles = secondaryFiles, formatOption = formatOption)
        }

      case womArray: WomArray =>
        womArray.value.toList.traverse(
          populateSecondaryFiles(_, secondaryFilesCoproduct, formatOption, parameterContext, expressionLib, ioFunctions)
        ).map(WomArray(_))

      case womValue: WomValue => womValue.valid
    }
  }

  /**
    * Checks if the file is compatible with a format.
    */
  def checkFormat(womMaybePopulatedFile: WomMaybePopulatedFile,
                  formatsOption: Option[List[String]],
                  schemaOption: Option[Schema]): ErrorOr[Unit] = {
    validate {
      for {
        schema <- schemaOption
        fileFormat <- womMaybePopulatedFile.formatOption
        formats <- formatsOption
      } yield {
        if (!formats.exists(schema.isSubClass(fileFormat, _)))
          throw new RuntimeException(s"$fileFormat is not compatible with ${formats.mkString(", ")}")
      }
      ()
    }
  }

  /**
    * Populates the contents if they aren't loaded already.
    */
  def maybeLoadContents(womMaybePopulatedFile: WomMaybePopulatedFile,
                        ioFunctionSet: IoFunctionSet,
                        loadContents: Boolean): ErrorOr[Option[String]] = {
    womMaybePopulatedFile.contentsOption match {
      case someContents@Some(_) => someContents.valid
      case None if !loadContents => None.valid
      case _ => FileParameter.load64KiB(womMaybePopulatedFile.value, ioFunctionSet).map(Option(_))
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
                     expressionLib: ExpressionLib,
                     ioFunctions: IoFunctionSet): ErrorOr[List[WomFile]] = {
    secondaryFilesOption
      .map(secondaryFiles(primaryWomFile, stringWomFileType, _, parameterContext, expressionLib, ioFunctions))
      .getOrElse(Nil.valid)
  }

  /**
    * Returns the list of secondary files for the primary file.
    */
  def secondaryFiles(primaryWomFile: WomFile,
                     stringWomFileType: WomFileType,
                     secondaryFiles: SecondaryFiles,
                     parameterContext: ParameterContext,
                     expressionLib: ExpressionLib,
                     ioFunctions: IoFunctionSet): ErrorOr[List[WomFile]] = {
    secondaryFiles
      .fold(SecondaryFilesPoly)
      .apply(primaryWomFile, stringWomFileType, parameterContext, expressionLib, ioFunctions)
  }

  type SecondaryFilesFunction = (WomFile, WomFileType, ParameterContext, ExpressionLib, IoFunctionSet) => ErrorOr[List[WomFile]]

  object SecondaryFilesPoly extends Poly1 {
    implicit def caseStringOrExpression: Case.Aux[StringOrExpression, SecondaryFilesFunction] = {
      at {
        _.fold(this)
      }
    }

    implicit def caseExpression: Case.Aux[Expression, SecondaryFilesFunction] = {
      at {
        expression =>
          (primaryWomFile, stringWomFileType, parameterContext, expressionLib, ioFunctions) =>
            File.secondaryExpressionFiles(primaryWomFile, stringWomFileType, expression, parameterContext, expressionLib, ioFunctions)
      }
    }

    implicit def caseString: Case.Aux[String, SecondaryFilesFunction] = {
      at {
        string =>
          (primaryWomFile, stringWomFileType, _, _, ioFunctions) =>
            File.secondaryStringFile(primaryWomFile, stringWomFileType, string, ioFunctions).map(List(_))
      }
    }

    implicit def caseArray: Case.Aux[Array[StringOrExpression], SecondaryFilesFunction] = {
      at {
        array =>
          (primaryWomFile, stringWomFileType, parameterContext, expressionLib, ioFunctions) =>
            val functions: List[SecondaryFilesFunction] = array.toList.map(_.fold(this))
            functions.flatTraverse(_ (primaryWomFile, stringWomFileType, parameterContext, expressionLib, ioFunctions))
      }
    }
  }

}
