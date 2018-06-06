package cwl

import cats.instances.future._
import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation._
import wom.expression.{FileEvaluation, IoFunctionSet}
import wom.expression.IoFunctionSet.{IoDirectory, IoFile}
import wom.types._
import wom.values._

import scala.concurrent.Await
import scala.concurrent.duration._

/** @see <a href="http://www.commonwl.org/v1.0/Workflow.html#CommandOutputBinding">CommandOutputBinding</a> */
case class CommandOutputBinding(glob: Option[Glob] = None,
                                loadContents: Option[Boolean] = None,
                                outputEval: Option[StringOrExpression] = None)

object CommandOutputBinding {
  /**
    * TODO: WOM: WOMFILE: Need to support returning globs for primary and secondary.
    *
    * Right now, when numerous glob files are created, the backends place each glob result into a different subdirectory.
    * Later, when the secondary files are being generated each of their globs are in different sub-directories.
    * Unfortunately, the secondary file resolution assumes that the secondaries are _next to_ the primaries, not in
    * sibling directories.
    *
    * While this is being worked on, instead of looking up globs return regular files until this can be fixed.
    */
  def isRegularFile(path: String): Boolean = {
    path match {
      case _ if path.contains("*") => false
      case _ if path.contains("?") => false
      case _ => true
    }
  }

  /**
    * Returns all the primary and secondary files that _will be_ created by this command output binding.
    */
  def getOutputWomFiles(inputValues: Map[String, WomValue],
                        outputWomType: WomType,
                        commandOutputBinding: CommandOutputBinding,
                        secondaryFilesOption: Option[SecondaryFiles],
                        ioFunctionSet: IoFunctionSet,
                        expressionLib: ExpressionLib): ErrorOr[Set[FileEvaluation]] = {
    val parameterContext = ParameterContext(ioFunctionSet, expressionLib, inputs = inputValues)

    /*
    CWL can output two types of "glob" path types:
    - File will be a glob path
    - Directory is a path that will be searched for files. It's assumed that it will not be a glob, but the spec is
      unclear and doesn't specifically say how Directories are listed. If the backend uses a POSIX
      `find <dir> -type f`, then a globbed path should still work.

    Either way create one of the two types as a return type that the backend should be looking for when the command
    completes.

    UPDATE:
    It turns out that secondary files can be directories. TODO: WOM: WOMFILE: Not sure what to do with that yet.
    https://github.com/common-workflow-language/common-workflow-language/blob/master/v1.0/v1.0/search.cwl#L42
    https://github.com/common-workflow-language/common-workflow-language/blob/master/v1.0/v1.0/index.py#L36-L37
     */
    def primaryPathsToWomFiles(primaryPaths: List[String]): ErrorOr[List[WomFile]] = {
      validate {
        primaryPaths map { primaryPath =>
          val outputWomFlatType = outputWomType match {
            case WomMaybeListedDirectoryType => WomUnlistedDirectoryType
            case WomArrayType(WomMaybeListedDirectoryType) => WomUnlistedDirectoryType
            case _ if isRegularFile(primaryPath) => WomSingleFileType
            case _ => WomGlobFileType
          }

          WomFile(outputWomFlatType, primaryPath)
        }
      }
    }

    def secondaryFilesToWomFiles(primaryWomFiles: List[WomFile], ioFunctionSet: IoFunctionSet): ErrorOr[List[WomFile]] = {
      primaryWomFiles.flatTraverse[ErrorOr, WomFile] { primaryWomFile =>
        FileParameter.secondaryFiles(primaryWomFile,
          primaryWomFile.womFileType, secondaryFilesOption, parameterContext, expressionLib, ioFunctionSet)
      }
    }

    for {
      primaryPaths <- GlobEvaluator.globs(commandOutputBinding.glob, parameterContext, expressionLib)
      primaryWomFiles <- primaryPathsToWomFiles(primaryPaths)
      // This sets optional = false arbitrarily for now as this code doesn't have the context to make that determination,
      // the caller can change this if necessary.
      primaryEvaluations = primaryWomFiles map { FileEvaluation(_, optional = false, secondary = false) }
      secondaryWomFiles <- secondaryFilesToWomFiles(primaryWomFiles, ioFunctionSet)
      secondaryEvaluations = secondaryWomFiles map { FileEvaluation(_, optional = false, secondary = true) }
    } yield (primaryEvaluations ++ secondaryEvaluations).toSet
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
    val parameterContext = ParameterContext(ioFunctionSet, expressionLib, inputs = inputValues)
    implicit val ec = ioFunctionSet.ec

    // 3. outputEval: pass in the primary files to an expression to generate our return value
    def evaluateWomValue(womFilesArray: WomArray): ErrorOr[WomValue] = {
      commandOutputBinding.outputEval match {
        case Some(StringOrExpression.String(string)) => WomString(string).valid
        case Some(StringOrExpression.Expression(expression)) =>
          val outputEvalParameterContext = parameterContext.copy(self = womFilesArray)
          ExpressionEvaluator.eval(expression, outputEvalParameterContext)
        case None =>
          womFilesArray.valid
      }
    }

    // Used to retrieve the file format to be injected into a file result.
    def formatOptionErrorOr = OutputParameter.format(formatCoproduct, parameterContext, expressionLib)

    // 4. secondaryFiles: just before returning the value, fill in the secondary files on the return value
    def populateSecondaryFiles(evaluatedWomValue: WomValue): ErrorOr[WomValue] = {
      for {
        formatOption <- formatOptionErrorOr
        womValue <- FileParameter.populateSecondaryFiles(
          evaluatedWomValue,
          secondaryFilesCoproduct,
          formatOption,
          parameterContext,
          expressionLib,
          ioFunctionSet
        )
      } yield womValue
    }

    // CWL tells us the type this output is expected to be. Attempt to coerce the actual output into this type.
    def coerceWomValue(populatedWomValue: WomValue): ErrorOr[WomValue] = {
      (outputWomType, populatedWomValue) match {

        case (womType: WomArrayType, womValue: WomArray) =>
          // Array -from- Array then coerce normally
          // NOTE: this evaluates to the same as the last case, but guards against the next cases accidentally matching
          womType.coerceRawValue(womValue).toErrorOr

        case (womType: WomArrayType, womValue) =>
          // Array -from- Single then coerce a single value to an array
          womType.coerceRawValue(womValue).toErrorOr.map(value => WomArray(womType, List(value)))

        case (womType, womValue: WomArray) if womValue.value.lengthCompare(1) == 0 =>
          // Single -from- Array then coerce the head (if there's only a head)
          womType.coerceRawValue(womValue.value.head).toErrorOr

        case (womType, womValue) =>
          // <other> to <other> then coerce normally
          womType.coerceRawValue(womValue).toErrorOr

      }
    }

    for {
      // 1. glob: get a list the globbed files as our primary files
      primaryPaths <- GlobEvaluator.globs(commandOutputBinding.glob, parameterContext, expressionLib)

      // 2. loadContents: load the contents of the primary files
      primaryAsDirectoryOrFiles <- primaryPaths.flatTraverse{
        loadPrimaryWithContents(ioFunctionSet, outputWomType, commandOutputBinding)
      }

      // Make globbed files absolute paths by prefixing them with the output dir if necessary
      absolutePaths = primaryAsDirectoryOrFiles.map(_.mapFile(ioFunctionSet.pathFunctions.relativeToHostCallRoot))
      
      // Load file size
      withFileSizes <- FileParameter.sync(absolutePaths.traverse(_.withSize(ioFunctionSet))).toErrorOr

      womFilesArray = WomArray(withFileSizes)

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
        loadDirectoryWithListing(ioFunctionSet, commandOutputBinding)(cwlPath).map(List(_))
      case WomMaybePopulatedFileType if isRegularFile(cwlPath) =>
        val globs = List(cwlPath)
        globs traverse loadFileWithContents(ioFunctionSet, commandOutputBinding)
      case WomMaybePopulatedFileType =>
        val globs: Seq[String] = Await.result(ioFunctionSet.glob(cwlPath), Duration.Inf).

        //TODO: HACK ALERT - DB: I am starting on ticket https://github.com/broadinstitute/cromwell/issues/3472 which will redeem me of this mortal sin.
          filterNot{s =>
            s.endsWith("rc.tmp") || s.endsWith("docker_cid") || s.endsWith("script") || s.endsWith("script.background") || s.endsWith("script.submit") || s.endsWith("stderr") || s.endsWith("stderr.background") || s.endsWith("stdout") || s.endsWith("stdout.background")
          }

        globs.toList traverse loadFileWithContents(ioFunctionSet, commandOutputBinding)
      case other => s"Program error: $other type was not expected".invalidNel
    }
  }

  /**
    * Loads a directory with the files listed, each file with contents populated.
    */
  private def loadDirectoryWithListing(ioFunctionSet: IoFunctionSet,
                                       commandOutputBinding: CommandOutputBinding)(path: String, visited: Vector[String] = Vector.empty): ErrorOr[WomMaybeListedDirectory] = {
    val listing = Await.result(ioFunctionSet.listDirectory(path)(visited), Duration.Inf).toList
    listing.traverse[ErrorOr, WomFile]({
      case IoFile(p) => loadFileWithContents(ioFunctionSet, commandOutputBinding)(p)
      case IoDirectory(p) => loadDirectoryWithListing(ioFunctionSet, commandOutputBinding)(p, visited :+ path)
    }).map({ listing =>
      WomMaybeListedDirectory(valueOption = Option(path), listingOption = Option(listing))
    })
  }

  /**
    * Loads a file at path reading 64KiB of data into the contents.
    */
  private def loadFileWithContents(ioFunctionSet: IoFunctionSet,
                                   commandOutputBinding: CommandOutputBinding)(path: String): ErrorOr[WomFile] = {
    val contentsOptionErrorOr = {
      if (commandOutputBinding.loadContents getOrElse false) {
        FileParameter.load64KiB(path, ioFunctionSet) map Option.apply
      } else {
        None.valid
      }
    }
    contentsOptionErrorOr map { contentsOption =>
      WomMaybePopulatedFile(valueOption = Option(path), contentsOption = contentsOption)
    }
  }
}
