package cwl

import cats.effect.IO
import cats.effect.IO._
import cats.instances.list._
import cats.syntax.functor._
import cats.syntax.traverse._
import cats.syntax.validated._
import cats.syntax.parallel._
import common.validation.ErrorOr._
import common.validation.IOChecked._
import common.validation.Validation._
import wom.expression.IoFunctionSet.{IoDirectory, IoFile}
import wom.expression.{FileEvaluation, IoFunctionSet}
import wom.types._
import wom.values._

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
                        expressionLib: ExpressionLib): IOChecked[Set[FileEvaluation]] = {
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

    def secondaryFilesToWomFiles(primaryWomFiles: List[WomFile], ioFunctionSet: IoFunctionSet): IOChecked[List[WomFile]] = {
      primaryWomFiles.flatTraverse[IOChecked, WomFile] { primaryWomFile =>
        FileParameter.secondaryFiles(primaryWomFile,
          primaryWomFile.womFileType, secondaryFilesOption, parameterContext, expressionLib, ioFunctionSet)
      }
    }

    for {
      primaryPaths <- GlobEvaluator.globs(commandOutputBinding.glob, parameterContext, expressionLib).toIOChecked
      primaryWomFiles <- primaryPathsToWomFiles(primaryPaths).toIOChecked
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
                             expressionLib: ExpressionLib): IOChecked[WomValue] = {
    import ioFunctionSet.cs

    val parameterContext = ParameterContext(ioFunctionSet, expressionLib, inputs = inputValues)

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
    def populateSecondaryFiles(evaluatedWomValue: WomValue): IOChecked[WomValue] = {
      for {
        formatOption <- formatOptionErrorOr.toIOChecked
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
      primaryPaths <- GlobEvaluator.globs(commandOutputBinding.glob, parameterContext, expressionLib).toIOChecked

      // 2. loadContents: load the contents of the primary files
      primaryAsDirectoryOrFiles <- primaryPaths.parTraverse[IOChecked, IOCheckedPar, List[WomFile]] {
        loadPrimaryWithContents(ioFunctionSet, outputWomType, commandOutputBinding)
      } map (_.flatten)

      // Make globbed files absolute paths by prefixing them with the output dir if necessary
      absolutePaths = primaryAsDirectoryOrFiles.map(_.mapFile(ioFunctionSet.pathFunctions.relativeToHostCallRoot))
      
      // Load file size
      withFileSizes <- absolutePaths.parTraverse[IOChecked, IOCheckedPar, WomFile](_.withSize(ioFunctionSet).to[IOChecked])

      womFilesArray = WomArray(withFileSizes)

      // 3. outputEval: pass in the primary files to an expression to generate our return value
      evaluatedWomValue <- evaluateWomValue(womFilesArray).toIOChecked

      // 4. secondaryFiles: just before returning the value, fill in the secondary files on the return value
      populatedWomValue <- populateSecondaryFiles(evaluatedWomValue)

      // CWL tells us the type this output is expected to be. Attempt to coerce the actual output into this type.
      coercedWomValue <- coerceWomValue(populatedWomValue).toIOChecked
    } yield coercedWomValue
  }

  /**
    * Given a cwl glob path and an output type, gets the listing using the ioFunctionSet, and optionally loads the
    * contents of the file(s).
    */
  private def loadPrimaryWithContents(ioFunctionSet: IoFunctionSet,
                                      outputWomType: WomType,
                                      commandOutputBinding: CommandOutputBinding)
                                     (cwlPath: String): IOChecked[List[WomFile]] = {
    import ioFunctionSet.cs

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
        loadFileWithContents(ioFunctionSet, commandOutputBinding)(cwlPath).to[IOChecked].map(List(_))
      case WomMaybePopulatedFileType =>
        //TODO: HACK ALERT - DB: I am starting on ticket https://github.com/broadinstitute/cromwell/issues/3092 which will redeem me of this mortal sin.
        val detritusFiles = List(
          "docker_cid",
          "gcs_delocalization.sh",
          "gcs_localization.sh",
          "gcs_transfer.sh",
          "rc.tmp",
          "script",
          "script.background",
          "script.submit",
          "stderr",
          "stderr.background",
          "stdout",
          "stdout.background",
        )
        val globs: IOChecked[Seq[String]] = 
          ioFunctionSet.glob(cwlPath).toIOChecked
              .map({
                _
                .filterNot{ s =>
                  detritusFiles exists s.endsWith
                }
              })

        globs.flatMap({ files =>
          files.toList.parTraverse[IOChecked, IOCheckedPar, WomFile](v => loadFileWithContents(ioFunctionSet, commandOutputBinding)(v).to[IOChecked])
        }) 
      case other => s"Program error: $other type was not expected".invalidIOChecked
    }
  }

  /**
    * Loads a directory with the files listed, each file with contents populated.
    */
  private def loadDirectoryWithListing(ioFunctionSet: IoFunctionSet,
                                       commandOutputBinding: CommandOutputBinding)(path: String, visited: Vector[String] = Vector.empty): IOChecked[WomMaybeListedDirectory] = {
    import ioFunctionSet.cs

    for {
      listing <- IO.fromFuture(IO { ioFunctionSet.listDirectory(path)(visited) }).to[IOChecked]
      loadedListing <- listing.toList.parTraverse[IOChecked, IOCheckedPar, WomFile]({
        case IoFile(p) => loadFileWithContents(ioFunctionSet, commandOutputBinding)(p).to[IOChecked]
        case IoDirectory(p) => loadDirectoryWithListing(ioFunctionSet, commandOutputBinding)(p, visited :+ path).widen
      })
    } yield WomMaybeListedDirectory(valueOption = Option(path), listingOption = Option(loadedListing))
  }

  /**
    * Loads a file at path reading 64KiB of data into the contents.
    */
  private def loadFileWithContents(ioFunctionSet: IoFunctionSet,
                                   commandOutputBinding: CommandOutputBinding)(path: String): IO[WomFile] = {
    val contentsOptionErrorOr = {
      if (commandOutputBinding.loadContents getOrElse false) {
        FileParameter.load64KiB(path, ioFunctionSet) map Option.apply
      } else {
        IO.pure(None)
      }
    }
    contentsOptionErrorOr map { contentsOption =>
      WomMaybePopulatedFile(valueOption = Option(path), contentsOption = contentsOption)
    }
  }
}
