package cwl

import cats.effect.IO
import cats.instances.list._
import cats.syntax.traverse._
import common.validation.ErrorOr._
import common.validation.IOChecked._
import common.validation.Validation._
import cwl.ontology.Schema
import shapeless.Poly1
import wom.expression.IoFunctionSet
import wom.types.{WomFileType, WomMaybePopulatedFileType}
import wom.values.{WomArray, WomFile, WomMaybePopulatedFile, WomValue}

object FileParameter {
  private val ReadLimit = Option(64 * 1024)

  def populateSecondaryFiles(womValue: WomValue,
                             secondaryFilesCoproduct: Option[SecondaryFiles],
                             formatOption: Option[String],
                             parameterContext: ParameterContext,
                             expressionLib: ExpressionLib,
                             ioFunctions: IoFunctionSet): IOChecked[WomValue] = {

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

      case womValue: WomValue => womValue.validIOChecked
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
                        loadContents: Boolean): IO[Option[String]] = {
    womMaybePopulatedFile.contentsOption match {
      case someContents@Some(_) => IO.pure(someContents)
      case None if !loadContents => IO.pure(None)
      case _ => FileParameter.load64KiB(womMaybePopulatedFile.value, ioFunctionSet).map(Option(_))
    }
  }

  def load64KiB(path: String, ioFunctionSet: IoFunctionSet): IO[String] = {
    IO.fromFuture(IO { ioFunctionSet.readFile(path, ReadLimit, failOnOverflow = false) })
  }

  /**
    * Returns the list of secondary files for the primary file.
    */
  def secondaryFiles(primaryWomFile: WomFile,
                     stringWomFileType: WomFileType,
                     secondaryFilesOption: Option[SecondaryFiles],
                     parameterContext: ParameterContext,
                     expressionLib: ExpressionLib,
                     ioFunctions: IoFunctionSet): IOChecked[List[WomFile]] = {
    secondaryFilesOption
      .map(secondaryFiles(primaryWomFile, stringWomFileType, _, parameterContext, expressionLib, ioFunctions))
      .getOrElse(List.empty[WomFile].validIOChecked)
  }

  /**
    * Returns the list of secondary files for the primary file.
    */
  def secondaryFiles(primaryWomFile: WomFile,
                     stringWomFileType: WomFileType,
                     secondaryFiles: SecondaryFiles,
                     parameterContext: ParameterContext,
                     expressionLib: ExpressionLib,
                     ioFunctions: IoFunctionSet): IOChecked[List[WomFile]] = {
    secondaryFiles
      .fold(SecondaryFilesPoly)
      .apply(primaryWomFile, stringWomFileType, parameterContext, expressionLib, ioFunctions)
  }

  type SecondaryFilesFunction = (WomFile, WomFileType, ParameterContext, ExpressionLib, IoFunctionSet) => IOChecked[List[WomFile]]

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
            File.secondaryExpressionFiles(primaryWomFile, stringWomFileType, expression, parameterContext, expressionLib, ioFunctions).toIOChecked
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
