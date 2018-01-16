package cwl

import cats.instances.list._
import cats.syntax.apply._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation._
import wom.expression.IoFunctionSet
import wom.types._
import wom.values._

import scala.concurrent.Await
import scala.concurrent.duration._

/** @see <a href="http://www.commonwl.org/v1.0/Workflow.html#CommandOutputBinding">CommandOutputBinding</a> */
case class CommandOutputBinding(glob: Option[Glob] = None,
                                loadContents: Option[Boolean] = None,
                                outputEval: Option[StringOrExpression] = None)

object CommandOutputBinding {
  val ReadLimit = Option(64 * 1024)

  /**
    * Returns all the primary and secondary files that _will be_ created by this command output binding.
    */
  def getOutputWomFiles(inputValues: Map[String, WomValue],
                        outputWomType: WomType,
                        commandOutputBinding: CommandOutputBinding,
                        secondaryFilesOption: Option[SecondaryFiles],
                        expressionLib: ExpressionLib): ErrorOr[Set[WomFile]] = {
    val parameterContext = ParameterContext(inputs = inputValues)

    /*
    CWL can output two types of "glob" path types:
    - File will be a glob path
    - Directory is a path that will be searched for files. It's assumed that it will not be a glob, but the spec is
      unclear and doesn't specifically say how Directories are listed. If the backend uses a POSIX
      `find <dir> -type f`, then a globbed path should still work.

    Either way create one of the two types as a return type that the backend should be looking for when the command
    completes.
     */
    val outputWomFlatType = outputWomType match {
      case WomMaybeListedDirectoryType => WomUnlistedDirectoryType
      case WomArrayType(WomMaybeListedDirectoryType) => WomUnlistedDirectoryType
      case _ => WomGlobFileType
    }

    for {
      primaryPaths <- GlobEvaluator.globs(commandOutputBinding.glob, parameterContext, expressionLib)
      primaryWomFiles <- outputWomFlatType match {
        case WomGlobFileType => primaryPaths.map(WomGlobFile).valid
        case WomUnlistedDirectoryType => primaryPaths.map(WomUnlistedDirectory).valid
        case other => s"Program error: $other type was not expected".invalidNel
      }
      secondaryWomFiles <- primaryWomFiles.flatTraverse[ErrorOr, WomFile] {
        CommandLineTool.CommandOutputParameter.secondaryFiles(_, secondaryFilesOption, parameterContext, expressionLib)
      }
    } yield (primaryWomFiles ++ secondaryWomFiles).toSet
  }

  /**
    * Generates an output wom value based on the specification in command output binding.
    *
    * Depending on the outputWomType, the following steps will be applied as specified in the CWL spec:
    * 1. glob: get a list the globbed files as our primary files
    * 2. loadContents: load the contents of the primary files
    * 3. outputEval: pass in the primary files to an expression to generate our return value
    * 4. secondaryFiles: just before returning the value, fill in the secondary files on the return value
    *
    * The result type will be coerced to the output type.
    */
  def generateOutputWomValue(inputValues: Map[String, WomValue],
                             ioFunctionSet: IoFunctionSet,
                             outputWomType: WomType,
                             commandOutputBinding: CommandOutputBinding,
                             secondaryFilesCoproduct: Option[SecondaryFiles],
                             formatCoproduct: Option[StringOrExpression],
                             expressionLib: ExpressionLib): ErrorOr[WomValue] = {
    val parameterContext = ParameterContext(inputs = inputValues)

    // 3. outputEval: pass in the primary files to an expression to generate our return value
    def evaluateWomValue(womFilesArray: WomArray): ErrorOr[WomValue] = {
      commandOutputBinding.outputEval match {
        case Some(StringOrExpression.String(string)) => WomString(string).valid
        case Some(StringOrExpression.Expression(expression)) =>
          val outputEvalParameterContext = parameterContext.copy(self = womFilesArray)
          ExpressionEvaluator.eval(expression, outputEvalParameterContext, expressionLib)
        case None =>
          womFilesArray.valid
      }
    }

    // Used to retrieve the file format to be injected into a file result.
    def formatOptionErrorOr = CommandLineTool.CommandOutputParameter.format(formatCoproduct, parameterContext, expressionLib)

    // 4. secondaryFiles: just before returning the value, fill in the secondary files on the return value
    def populateSecondaryFiles(evaluatedWomValue: WomValue): ErrorOr[WomValue] = {
      evaluatedWomValue match {
        case womMaybePopulatedFile: WomMaybePopulatedFile =>
          val secondaryFilesErrorOr = CommandLineTool.CommandOutputParameter.secondaryFiles(
            womMaybePopulatedFile,
            secondaryFilesCoproduct,
            parameterContext,
            expressionLib
          )

          (secondaryFilesErrorOr, formatOptionErrorOr) mapN { (secondaryFiles, formatOption) =>
            womMaybePopulatedFile.copy(secondaryFiles = secondaryFiles, formatOption = formatOption)
          }

        case womArray: WomArray if womArray.womType.memberType == WomMaybePopulatedFileType =>
          formatOptionErrorOr map { formatOption =>
            womArray.map(_.asInstanceOf[WomMaybePopulatedFile].copy(formatOption = formatOption))
          }

        case womValue: WomValue => womValue.valid
      }
    }

    // CWL tells us the type this output is expected to be. Attempt to coerce the actual output into this type.
    def coerceWomValue(populatedWomValue: WomValue): ErrorOr[WomValue] = {
      (outputWomType, populatedWomValue) match {
        case (womType: WomArrayType, womValue: WomArray) => womType.coerceRawValue(womValue).toErrorOr
        case (womType: WomArrayType, womValue) =>
          // Coerce a single value to an array
          womType.coerceRawValue(womValue).toErrorOr.map(value => WomArray(womType, List(value)))
        case (womType, womValue: WomArray) if womValue.value.lengthCompare(1) == 0 =>
          // Coerce an array to a single value
          womType.coerceRawValue(womValue.value.head).toErrorOr
        case (womType, womValue) => womType.coerceRawValue(womValue).toErrorOr
      }
    }

    for {
      // 1. glob: get a list the globbed files as our primary files
      primaryPaths <- GlobEvaluator.globs(commandOutputBinding.glob, parameterContext, expressionLib)

      // 2. loadContents: load the contents of the primary files
      primaryAsDirectoryOrFiles <- primaryPaths.flatTraverse[ErrorOr, WomFile] {
        loadPrimaryWithContents(ioFunctionSet, outputWomType, commandOutputBinding)
      }
      womFilesArray = WomArray(primaryAsDirectoryOrFiles)

      // 3. outputEval: pass in the primary files to an expression to generate our return value
      evaluatedWomValue <- evaluateWomValue(womFilesArray)

      // 4. secondaryFiles: just before returning the value, fill in the secondary files on the return value
      populatedWomValue <- populateSecondaryFiles(evaluatedWomValue)

      // CWL tells us the type this output is expected to be. Attempt to coerce the actual output into this type.
      coercedWomValue <- coerceWomValue(populatedWomValue)
    } yield coercedWomValue
  }

  /**
    * Given a cwl glob path and an output type, gets the listing using the ioFunctionSet, and optionally loads the
    * contents of the file(s).
    */
  private def loadPrimaryWithContents(ioFunctionSet: IoFunctionSet,
                                      outputWomType: WomType,
                                      commandOutputBinding: CommandOutputBinding)
                                     (cwlPath: String): ErrorOr[List[WomFile]] = {
    /*
    For each file matched in glob, read up to the first 64 KiB of text from the file and place it in the contents field
    of the file object for manipulation by outputEval.
     */
    val womMaybeListedDirectoryOrFileType = outputWomType match {
      case WomMaybeListedDirectoryType => WomMaybeListedDirectoryType
      case WomArrayType(WomMaybeListedDirectoryType) => WomMaybeListedDirectoryType
      case _ => WomMaybePopulatedFileType
    }
    womMaybeListedDirectoryOrFileType match {
      case WomMaybeListedDirectoryType =>
        // Even if multiple directories are _somehow_ requested, a single flattened directory is returned.
        loadDirectoryWithListing(cwlPath, ioFunctionSet, commandOutputBinding).map(List(_))
      case WomMaybePopulatedFileType =>
        val globs = Await.result(ioFunctionSet.glob(cwlPath), Duration.Inf)
        globs.toList traverse loadFileWithContents(ioFunctionSet, commandOutputBinding)
      case other => s"Program error: $other type was not expected".invalidNel
    }
  }

  /**
    * Loads a directory with the files listed, each file with contents populated.
    */
  private def loadDirectoryWithListing(path: String,
                                       ioFunctionSet: IoFunctionSet,
                                       commandOutputBinding: CommandOutputBinding): ErrorOr[WomMaybeListedDirectory] = {
    val listing = Await.result(ioFunctionSet.listAllFilesUnderDirectory(path), Duration.Inf)
    listing.toList traverse loadFileWithContents(ioFunctionSet, commandOutputBinding) map { listing =>
      WomMaybeListedDirectory(valueOption = Option(path), listingOption = Option(listing))
    }
  }

  /**
    * Loads a file at path reading 64KiB of data into the contents.
    */
  private def loadFileWithContents(ioFunctionSet: IoFunctionSet,
                                   commandOutputBinding: CommandOutputBinding)(path: String): ErrorOr[WomFile] = {
    val contentsOptionErrorOr = {
      if (commandOutputBinding.loadContents getOrElse false) {
        load64KiB(path, ioFunctionSet) map Option.apply
      } else {
        None.valid
      }
    }
    contentsOptionErrorOr map { contentsOption =>
      WomMaybePopulatedFile(valueOption = Option(path), contentsOption = contentsOption)
    }
  }

  private def load64KiB(path: String, ioFunctionSet: IoFunctionSet): ErrorOr[String] = {
    // TODO: WOM: propagate Future (or IO?) signature
    validate {
      val content = ioFunctionSet.readFile(path, ReadLimit, failOnOverflow = false)
      Await.result(content, 60.seconds)
    }
  }
}
