package cromwell.backend.standard

import java.io.IOException

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.event.LoggingReceive
import cats.instances.list._
import cats.instances.option._
import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.functor._
import cats.syntax.traverse._
import common.Checked
import common.exception.MessageAggregation
import common.util.StringUtil._
import common.util.TryUtil
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import common.validation.IOChecked._
import common.validation.Validation._
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, JobAbortedResponse, JobReconnectionNotSupportedException}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.OutputEvaluator._
import cromwell.backend.SlowJobWarning.{WarnAboutSlownessAfter, WarnAboutSlownessIfNecessary}
import cromwell.backend._
import cromwell.backend.async.AsyncBackendJobExecutionActor._
import cromwell.backend.async._
import cromwell.backend.standard.StandardAdHocValue._
import cromwell.backend.validation._
import cromwell.core.io.{AsyncIoActorClient, DefaultIoCommandBuilder, IoCommandBuilder}
import cromwell.core.path.Path
import cromwell.core.{CromwellAggregatedException, CromwellFatalExceptionMarker, ExecutionEvent, StandardPaths}
import cromwell.services.keyvalue.KeyValueServiceActor._
import cromwell.services.keyvalue.KvClient
import cromwell.services.metadata.CallMetadataKeys
import mouse.all._
import net.ceedubs.ficus.Ficus._
import org.apache.commons.lang3.StringUtils
import shapeless.Coproduct
import wom.callable.{AdHocValue, CommandTaskDefinition, ContainerizedInputExpression, RuntimeEnvironment}
import wom.expression.WomExpression
import wom.graph.LocalName
import wom.values._
import wom.{CommandSetupSideEffectFile, InstantiatedCommand, WomFileMapper}

import scala.concurrent._
import scala.concurrent.duration._
import scala.util.{Failure, Success, Try}

trait StandardAsyncExecutionActorParams extends StandardJobExecutionActorParams {
  /** The promise that will be completed when the async run is complete. */
  def completionPromise: Promise[BackendJobExecutionResponse]
}

case class DefaultStandardAsyncExecutionActorParams
(
  override val jobIdKey: String,
  override val serviceRegistryActor: ActorRef,
  override val ioActor: ActorRef,
  override val jobDescriptor: BackendJobDescriptor,
  override val configurationDescriptor: BackendConfigurationDescriptor,
  override val backendInitializationDataOption: Option[BackendInitializationData],
  override val backendSingletonActorOption: Option[ActorRef],
  override val completionPromise: Promise[BackendJobExecutionResponse],
  override val minimumRuntimeSettings: MinimumRuntimeSettings
) extends StandardAsyncExecutionActorParams

/**
  * An extension of the generic AsyncBackendJobExecutionActor providing a standard abstract implementation of an
  * asynchronous polling backend.
  *
  * Backends supported by this trait will all have common behavior. If a backend implementor wishes to provide a custom
  * behavior, one should instead implement the various methods in `AsyncBackendJobExecutionActor`.
  *
  * NOTE: Unlike the parent trait `AsyncBackendJobExecutionActor`, this trait is subject to even more frequent updates
  * as the common behavior among the backends adjusts in unison.
  */
trait StandardAsyncExecutionActor
  extends AsyncBackendJobExecutionActor with StandardCachingActorHelper with AsyncIoActorClient with KvClient with SlowJobWarning {
  this: Actor with ActorLogging with BackendJobLifecycleActor =>

  override lazy val ioCommandBuilder: IoCommandBuilder = DefaultIoCommandBuilder

  val SIGTERM = 143
  val SIGINT = 130
  val SIGKILL = 137

  /** The type of the run info when a job is started. */
  type StandardAsyncRunInfo

  /** The type of the run status returned during each poll. */
  type StandardAsyncRunState

  /** Should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be carried around
    * in the state type.
    */
  def statusEquivalentTo(thiz: StandardAsyncRunState)(that: StandardAsyncRunState): Boolean

  /** The pending execution handle for each poll. */
  type StandardAsyncPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, StandardAsyncRunState]

  /** Standard set of parameters passed to the backend. */
  def standardParams: StandardAsyncExecutionActorParams

  override lazy val configurationDescriptor: BackendConfigurationDescriptor = standardParams.configurationDescriptor

  override lazy val completionPromise: Promise[BackendJobExecutionResponse] = standardParams.completionPromise

  override lazy val ioActor = standardParams.ioActor

  /** Backend initialization data created by the factory initializer. */
  override lazy val backendInitializationDataOption: Option[BackendInitializationData] =
    standardParams.backendInitializationDataOption

  /** @see [[StandardAsyncExecutionActorParams.serviceRegistryActor]] */
  override lazy val serviceRegistryActor: ActorRef = standardParams.serviceRegistryActor

  /** @see [[StandardAsyncExecutionActorParams.jobDescriptor]] */
  override lazy val jobDescriptor: BackendJobDescriptor = standardParams.jobDescriptor

  /** @see [[StandardAsyncExecutionActorParams.jobIdKey]] */
  lazy val jobIdKey: String = standardParams.jobIdKey

  /** @see [[Command.instantiate]] */
  lazy val backendEngineFunctions: StandardExpressionFunctions =
    standardInitializationData.expressionFunctions(jobPaths, standardParams.ioActor, ec)

  lazy val scriptEpilogue = configurationDescriptor.backendConfig.as[Option[String]]("script-epilogue").getOrElse("sync")

  lazy val temporaryDirectory = configurationDescriptor.backendConfig.getOrElse(
    path = "temporary-directory",
    default = s"""$$(mkdir -p "${runtimeEnvironment.tempPath}" && echo "${runtimeEnvironment.tempPath}")"""
  )

  def preProcessWomFile(womFile: WomFile): WomFile = womFile

  /** @see [[Command.instantiate]] */
  final def commandLinePreProcessor(inputs: WomEvaluatedCallInputs): Try[WomEvaluatedCallInputs] = {
    val map = inputs map { case (k, v) => k -> WomFileMapper.mapWomFiles(preProcessWomFile, inputsToNotLocalize)(v) }
    TryUtil.sequenceMap(map).
      recoverWith {
        case e => Failure(new IOException(e.getMessage) with CromwellFatalExceptionMarker)
      }
  }

  final lazy val localizedInputs: Try[WomEvaluatedCallInputs] = commandLinePreProcessor(jobDescriptor.evaluatedTaskInputs)

  /**
    * Maps WomFile to a local path, for use in the commandLineValueMapper.
    */
  def mapCommandLineWomFile(womFile: WomFile): WomFile =
    womFile.mapFile(workflowPaths.buildPath(_).pathAsString)

  // This is a trickery to allow the mapCommandLineWomFile above to call isAdHoc file without creating an
  // infinite recursion. This should really go away when we finally can have a sane implementation
  // keeping track of the paths cleanly without so many value mappers
  def mapCommandLineJobInputWomFile(womFile: WomFile): WomFile = mapCommandLineWomFile(womFile)

  // Allows backends to signal to the StandardAsyncExecutionActor that there's a set of input files which
  // they don't intend to localize for this task.
  def inputsToNotLocalize: Set[WomFile] = Set.empty

  /** @see [[Command.instantiate]] */
  final lazy val commandLineValueMapper: WomValue => WomValue = {
    womValue => WomFileMapper.mapWomFiles(mapCommandLineWomFile, inputsToNotLocalize)(womValue).get
  }

  /** @see [[Command.instantiate]] */
  final lazy val commandLineJobInputValueMapper: WomValue => WomValue = {
    womValue => WomFileMapper.mapWomFiles(mapCommandLineJobInputWomFile, inputsToNotLocalize)(womValue).get
  }

  lazy val jobShell: String = configurationDescriptor.backendConfig.getOrElse("job-shell",
    configurationDescriptor.globalConfig.getOrElse("system.job-shell", "/bin/bash"))

  lazy val abbreviateCommandLength: Int = configurationDescriptor.backendConfig.getOrElse("abbreviate-command-length",
    configurationDescriptor.globalConfig.getOrElse("system.abbreviate-command-length", 0))

  /**
    * The local path where the command will run.
    */
  lazy val commandDirectory: Path = jobPaths.callExecutionRoot

  /**
    * Returns the shell scripting for finding all files listed within a directory.
    *
    * @param directoryFiles The directories.
    * @return The shell scripting.
    */
  def directoryScripts(directoryFiles: Traversable[WomUnlistedDirectory]): String =
    directoryFiles map directoryScript mkString "\n"

  /**
    * Returns the shell scripting for finding all files listed within a directory.
    *
    * @param unlistedDirectory The directory.
    * @return The shell scripting.
    */
  def directoryScript(unlistedDirectory: WomUnlistedDirectory): String = {
    val directoryPath = unlistedDirectory.value.ensureUnslashed
    val directoryList = directoryPath + ".list"

    s"""|# list all the files that match the directory into a file called directory.list
        |find $directoryPath -type f > $directoryList
        |""".stripMargin
  }

  /**
    * The local parent directory of the glob file. By default this is the same as the commandDirectory.
    *
    * In some cases, to allow the hard linking by ln to operate, a different mount point must be returned.
    *
    * @param wdlGlobFile The glob.
    * @return The parent directory for writing the wdl glob.
    */
  def globParentDirectory(wdlGlobFile: WomGlobFile): Path = commandDirectory

  /**
    * Returns the shell scripting for hard linking the glob results using ln.
    *
    * @param globFiles The globs.
    * @return The shell scripting.
    */
  def globScripts(globFiles: Traversable[WomGlobFile]): String =
    globFiles map globScript mkString "\n"

  /**
    * Returns the shell scripting for linking a glob results file.
    *
    * @param globFile The glob.
    * @return The shell scripting.
    */
  def globScript(globFile: WomGlobFile): String = {
    val parentDirectory = globParentDirectory(globFile)
    val globDir = GlobFunctions.globName(globFile.value)
    val globDirectory = parentDirectory./(globDir)
    val globList = parentDirectory./(s"$globDir.list")
    val controlFileName = "cromwell_glob_control_file"
    val absoluteGlobValue = commandDirectory.resolve(globFile.value).pathAsString
    val globLinkCommand: String = configurationDescriptor.backendConfig.getAs[String]("glob-link-command")
      .map("( " + _ + " )")
      .getOrElse("( ln -L GLOB_PATTERN GLOB_DIRECTORY 2> /dev/null ) || ( ln GLOB_PATTERN GLOB_DIRECTORY )")
      .replaceAll("GLOB_PATTERN", absoluteGlobValue)
      .replaceAll("GLOB_DIRECTORY", globDirectory.pathAsString)

    val controlFileContent =
      """This file is used by Cromwell to allow for globs that would not match any file.
        |By its presence it works around the limitation of some backends that do not allow empty globs.
        |Regardless of the outcome of the glob, this file will not be part of the final list of globbed files.
      """.stripMargin

    s"""|# make the directory which will keep the matching files
        |mkdir $globDirectory
        |
        |# create the glob control file that will allow for the globbing to succeed even if there is 0 match
        |echo "${controlFileContent.trim}" > $globDirectory/$controlFileName
        |
        |# hardlink or symlink all the files into the glob directory
        |$globLinkCommand
        |
        |# list all the files (except the control file) that match the glob into a file called glob-[md5 of glob].list
        |ls -1 $globDirectory | grep -v $controlFileName > $globList
        |""".stripMargin
  }

  /** Any custom code that should be run within commandScriptContents before the instantiated command. */
  def scriptPreamble: String = ""

  def cwd: Path = commandDirectory
  def rcPath: Path = cwd./(jobPaths.returnCodeFilename)

  // The standard input filename can be as ephemeral as the execution: the name needs to match the expectations of
  // the command, but the standard input file will never be accessed after the command completes. standard output and
  // error on the other hand will be accessed and the names of those files need to be known to be delocalized and read
  // by the backend.
  //
  // All of these redirections can be expressions and will be given "value-mapped" versions of their inputs, which
  // means that if any fully qualified filenames are generated they will be container paths. As mentioned above this
  // doesn't matter for standard input since that's ephemeral to the container, but standard output and error paths
  // will need to be mapped back to host paths for use outside the command script.
  //
  // Absolutize any redirect and overridden paths. All of these files must have absolute paths since the command script
  // references them outside a (cd "execution dir"; ...) subshell. The default names are known to be relative paths,
  // the names from redirections may or may not be relative.
  private def absolutizeContainerPath(path: String): String = {
    if (path.startsWith(cwd.pathAsString)) path else cwd.resolve(path).pathAsString
  }

  def executionStdin: Option[String] = instantiatedCommand.evaluatedStdinRedirection map absolutizeContainerPath
  def executionStdout: String = instantiatedCommand.evaluatedStdoutOverride.getOrElse(jobPaths.defaultStdoutFilename) |> absolutizeContainerPath
  def executionStderr: String = instantiatedCommand.evaluatedStderrOverride.getOrElse(jobPaths.defaultStderrFilename) |> absolutizeContainerPath

  /*
   * Ensures the standard paths are correct w.r.t overridden paths. This is called in two places: when generating the command and
   * before sending the results back to the engine. The latter is to cover for the case of a restart where the command is not re-generated,
   * but the standard paths are still being used in the Success response to the engine. In that response paths need to be correct, hence
   * the call to this method. The var is to avoid doing the work twice: in non restarted runs the command is generated and so there is no need
   * to re-do this before sending the response.
   */
  private var jobPathsUpdated: Boolean = false
  private def updateJobPaths() = if (!jobPathsUpdated) {
    // .get's are safe on stdout and stderr after falling back to default names above.
    jobPaths.standardPaths = StandardPaths(
      output = hostPathFromContainerPath(executionStdout),
      error = hostPathFromContainerPath(executionStderr)
    )
    // Re-publish stdout and stderr paths that were possibly just updated.
    tellMetadata(jobPaths.standardOutputAndErrorPaths)
    jobPathsUpdated = true
  }

  def hostPathFromContainerPath(string: String): Path = {
    val cwdString = cwd.pathAsString.ensureSlashed
    val relativePath = string.stripPrefix(cwdString)
    jobPaths.callExecutionRoot.resolve(relativePath)
  }

  /** A bash script containing the custom preamble, the instantiated command, and output globbing behavior. */
  def commandScriptContents: ErrorOr[String] = {
    val commandString = instantiatedCommand.commandString
    val commandStringAbbreviated = StringUtils.abbreviateMiddle(commandString, "...", abbreviateCommandLength)
    jobLogger.info(s"`$commandStringAbbreviated`")
    tellMetadata(Map(CallMetadataKeys.CommandLine -> commandStringAbbreviated))

    val cwd = commandDirectory
    val rcPath = cwd./(jobPaths.returnCodeFilename)
    updateJobPaths()

    // The standard input redirection gets a '<' that standard output and error do not since standard input is only
    // redirected if requested. Standard output and error are always redirected but not necessarily to the default
    // stdout/stderr file names.
    val stdinRedirection = executionStdin.map("< " + _.shellQuote).getOrElse("")
    val stdoutRedirection = executionStdout.shellQuote
    val stderrRedirection = executionStderr.shellQuote
    val rcTmpPath = rcPath.plusExt("tmp")

    val errorOrDirectoryOutputs: ErrorOr[List[WomUnlistedDirectory]] =
      backendEngineFunctions.findDirectoryOutputs(call, jobDescriptor)

    val errorOrGlobFiles: ErrorOr[List[WomGlobFile]] =
      backendEngineFunctions.findGlobOutputs(call, jobDescriptor)

    lazy val environmentVariables = instantiatedCommand.environmentVariables map { case (k, v) => s"""export $k="$v"""" } mkString("", "\n", "\n")

    val home = jobDescriptor.taskCall.callable.homeOverride.map { _ (runtimeEnvironment) }.getOrElse("$HOME")
    val shortId = jobDescriptor.workflowDescriptor.id.shortString
    // Give the out and error FIFO variables names that are unlikely to conflict with anything the user is doing.
    val (out, err) = (s"out$shortId", s"err$shortId")

    val dockerOutputDir = jobDescriptor.taskCall.callable.dockerOutputDirectory map { d =>
      s"ln -s $cwd $d"
    } getOrElse ""

    // Only adjust the temporary directory permissions if this is executing under Docker.
    val tmpDirPermissionsAdjustment = if (isDockerRun) s"""chmod 777 "$$tmpDir"""" else ""

    val emptyDirectoryFillCommand: String = configurationDescriptor.backendConfig.getAs[String]("empty-dir-fill-command")
      .getOrElse(
        s"""(
           |# add a .file in every empty directory to facilitate directory delocalization on the cloud
           |cd $cwd
           |find . -type d -exec sh -c '[ -z "$$(ls -A '"'"'{}'"'"')" ] && touch '"'"'{}'"'"'/.file' \\;
           |)""".stripMargin)

    // The `tee` trickery below is to be able to redirect to known filenames for CWL while also streaming
    // stdout and stderr for PAPI to periodically upload to cloud storage.
    // https://stackoverflow.com/questions/692000/how-do-i-write-stderr-to-a-file-while-using-tee-with-a-pipe
    (errorOrDirectoryOutputs, errorOrGlobFiles).mapN((directoryOutputs, globFiles) =>
    s"""|#!$jobShell
        |DOCKER_OUTPUT_DIR_LINK
        |cd ${cwd.pathAsString}
        |tmpDir=$temporaryDirectory
        |$tmpDirPermissionsAdjustment
        |export _JAVA_OPTIONS=-Djava.io.tmpdir="$$tmpDir"
        |export TMPDIR="$$tmpDir"
        |export HOME="$home"
        |(
        |cd ${cwd.pathAsString}
        |SCRIPT_PREAMBLE
        |)
        |$out="$${tmpDir}/out.$$$$" $err="$${tmpDir}/err.$$$$"
        |mkfifo "$$$out" "$$$err"
        |trap 'rm "$$$out" "$$$err"' EXIT
        |tee $stdoutRedirection < "$$$out" &
        |tee $stderrRedirection < "$$$err" >&2 &
        |(
        |cd ${cwd.pathAsString}
        |ENVIRONMENT_VARIABLES
        |INSTANTIATED_COMMAND
        |) $stdinRedirection > "$$$out" 2> "$$$err"
        |echo $$? > $rcTmpPath
        |$emptyDirectoryFillCommand
        |(
        |cd ${cwd.pathAsString}
        |SCRIPT_EPILOGUE
        |${globScripts(globFiles)}
        |${directoryScripts(directoryOutputs)}
        |)
        |mv $rcTmpPath $rcPath
        |""".stripMargin
      .replace("SCRIPT_PREAMBLE", scriptPreamble)
      .replace("ENVIRONMENT_VARIABLES", environmentVariables)
      .replace("INSTANTIATED_COMMAND", commandString)
      .replace("SCRIPT_EPILOGUE", scriptEpilogue)
      .replace("DOCKER_OUTPUT_DIR_LINK", dockerOutputDir))
  }

  def runtimeEnvironmentPathMapper(env: RuntimeEnvironment): RuntimeEnvironment = {
    def localize(path: String): String = (WomSingleFile(path) |> commandLineValueMapper).valueString
    env.copy(outputPath = env.outputPath |> localize, tempPath = env.tempPath |> localize)
  }

  lazy val runtimeEnvironment = {
    RuntimeEnvironmentBuilder(jobDescriptor.runtimeAttributes, jobPaths)(standardParams.minimumRuntimeSettings) |> runtimeEnvironmentPathMapper
  }

  /**
    * By default, ad hoc values get localized to the call directory.
    * This way if running locally with docker they get mounted with the rest of the inputs in the container.
    * The PAPI backend overrides this to a noop since the localization happens on the VM directly, so there's no need
    * for this extra localization step.
    *
    * Maybe this should be the other way around: the default implementation is noop and SFS / TES override it ?
    */
  lazy val localizeAdHocValues: List[AdHocValue] => ErrorOr[List[StandardAdHocValue]] = { adHocValues =>
    import cats.instances.future._

    // Localize an adhoc file to the callExecutionRoot as needed
    val localize: (AdHocValue, Path) => Future[LocalizedAdHocValue] = { (adHocValue, file) =>
      val actualName = adHocValue.alternativeName.getOrElse(file.name)
      val finalPath = jobPaths.callExecutionRoot / actualName
      // First check that it's not already there under execution root
      asyncIo.existsAsync(finalPath) flatMap {
        // If it's not then copy it
        case false => asyncIo.copyAsync(file, finalPath) as { LocalizedAdHocValue(adHocValue, finalPath) }
        case true => Future.successful(LocalizedAdHocValue(adHocValue, finalPath))
      }
    }

    adHocValues.traverse[ErrorOr, (AdHocValue, Path)]({ adHocValue =>
      // Build an actionable Path from the ad hoc file
      getPath(adHocValue.womValue.value).toErrorOr.map(adHocValue -> _)
    })
      // Localize the values if necessary
      .map(_.traverse[Future, LocalizedAdHocValue](localize.tupled)).toEither
      // TODO: Asynchronify
      // This is obviously sad but turning it into a Future has earth-shattering consequences, so synchronizing it for now
      .flatMap(future => Try(Await.result(future, 1.hour)).toChecked)
      .map(_.map(Coproduct[StandardAdHocValue](_)))
      .toValidated
  }

  def adHocValueToCommandSetupSideEffectFile(adHocValue: StandardAdHocValue) = adHocValue match {
    case AsAdHocValue(AdHocValue(womValue, alternativeName, _)) =>
      CommandSetupSideEffectFile(womValue, alternativeName)
    case AsLocalizedAdHocValue(LocalizedAdHocValue(AdHocValue(womValue, alternativeName, _), _)) =>
      CommandSetupSideEffectFile(womValue, alternativeName)
  }

  lazy val evaluatedAdHocFiles: ErrorOr[List[AdHocValue]] = {
    val callable = jobDescriptor.taskCall.callable

    /*
      * Try to map the command line values.
      *
      * May not work as the commandLineValueMapper was originally meant to modify paths in the later stages of command
      * line instantiation. However, due to [[AdHocValue]] support the command line instantiation itself currently needs
      * to use the commandLineValue mapper. So the commandLineValueMapper is attempted first, and if that fails then
      * returns the original womValue.
      */
    def tryCommandLineValueMapper(womValue: WomValue): WomValue = {
      Try(commandLineJobInputValueMapper(womValue)).getOrElse(womValue)
    }

    val unmappedInputs: Map[String, WomValue] = jobDescriptor.evaluatedTaskInputs.map({
      case (inputDefinition, womValue) => inputDefinition.localName.value -> womValue
    })

    val mappedInputs: Checked[Map[String, WomValue]] = localizedInputs.toErrorOr.map(
      _.map({
        case (inputDefinition, value) => inputDefinition.localName.value -> tryCommandLineValueMapper(value)
      })
    ).toEither

    val evaluateAndInitialize = (containerizedInputExpression: ContainerizedInputExpression) => for {
      mapped <- mappedInputs
      evaluated <- containerizedInputExpression.evaluate(unmappedInputs, mapped, backendEngineFunctions).toChecked
      initialized <- evaluated.traverse[IOChecked, AdHocValue]({ adHocValue =>
        adHocValue.womValue.initialize(backendEngineFunctions).map({
          case file: WomFile => adHocValue.copy(womValue = file)
          case _ => adHocValue
        })
      }).toChecked
    } yield initialized

    callable.adHocFileCreation.toList
      .flatTraverse[ErrorOr, AdHocValue](evaluateAndInitialize.andThen(_.toValidated))
  }

  lazy val localizedAdHocValues: ErrorOr[List[StandardAdHocValue]] = evaluatedAdHocFiles.toEither
    .flatMap(localizeAdHocValues.andThen(_.toEither))
    .toValidated

  protected def asAdHocFile(womFile: WomFile) = evaluatedAdHocFiles map { _.find({
    case AdHocValue(file, _, _) => file.value == womFile.value
  })
  } getOrElse None

  protected def isAdHocFile(womFile: WomFile) = asAdHocFile(womFile).isDefined

  /** The instantiated command. */
  lazy val instantiatedCommand: InstantiatedCommand = {
    val callable = jobDescriptor.taskCall.callable

    // Replace input files with the ad hoc updated version
    def adHocFilePreProcessor(in: WomEvaluatedCallInputs): Try[WomEvaluatedCallInputs] = {
      localizedAdHocValues.toTry("Error evaluating ad hoc files") map { adHocFiles =>
        in map {
          case (inputDefinition, originalWomValue) =>
            inputDefinition -> adHocFiles.collectFirst({
              case AsAdHocValue(AdHocValue(originalWomFile, _, Some(inputName))) if inputName == inputDefinition.localName.value => originalWomFile
              case AsLocalizedAdHocValue(LocalizedAdHocValue(AdHocValue(originalWomFile, _, Some(inputName)), localizedPath)) if inputName == inputDefinition.localName.value =>
                originalWomFile.mapFile(_ => localizedPath.pathAsString)
            }).getOrElse(originalWomValue)
        }
      }
    }

    // Gets the inputs that will be mutated by instantiating the command.
    def mutatingPreProcessor(in: WomEvaluatedCallInputs): Try[WomEvaluatedCallInputs] = {
      for {
        commandLineProcessed <- localizedInputs
        adHocProcessed <- adHocFilePreProcessor(commandLineProcessed)
      } yield adHocProcessed
    }

    val instantiatedCommandValidation = Command.instantiate(
      jobDescriptor,
      backendEngineFunctions,
      mutatingPreProcessor,
      commandLineValueMapper,
      runtimeEnvironment
    )

    def makeStringKeyedMap(list: List[(LocalName, WomValue)]): Map[String, WomValue] = list.toMap map { case (k, v) => k.value -> v }

    val command = instantiatedCommandValidation flatMap { instantiatedCommand =>
      val valueMappedPreprocessedInputs = instantiatedCommand.valueMappedPreprocessedInputs |> makeStringKeyedMap

      val adHocFileCreationSideEffectFiles: ErrorOr[List[CommandSetupSideEffectFile]] = localizedAdHocValues map { _.map(adHocValueToCommandSetupSideEffectFile) }

      def evaluateEnvironmentExpression(nameAndExpression: (String, WomExpression)): ErrorOr[(String, String)] = {
        val (name, expression) = nameAndExpression
        expression.evaluateValue(valueMappedPreprocessedInputs, backendEngineFunctions) map { name -> _.valueString }
      }

      val environmentVariables = callable.environmentExpressions.toList traverse evaluateEnvironmentExpression

      // Build a list of functions from a CommandTaskDefinition to an Option[WomExpression] representing a possible
      // redirection or override of the filename of a redirection. Evaluate that expression if present and stringify.
      val List(stdinRedirect, stdoutOverride, stderrOverride) = List[CommandTaskDefinition => Option[WomExpression]](
        _.stdinRedirection, _.stdoutOverride, _.stderrOverride) map {
        _.apply(callable).traverse[ErrorOr, String] { _.evaluateValue(valueMappedPreprocessedInputs, backendEngineFunctions) map { _.valueString} }
      }

      (adHocFileCreationSideEffectFiles, environmentVariables, stdinRedirect, stdoutOverride, stderrOverride) mapN {
        (adHocFiles, env, in, out, err) =>
          instantiatedCommand.copy(
            createdFiles = instantiatedCommand.createdFiles ++ adHocFiles, environmentVariables = env.toMap,
            evaluatedStdinRedirection = in, evaluatedStdoutOverride = out, evaluatedStderrOverride = err)
      }: ErrorOr[InstantiatedCommand]
    }

    // TODO CWL: Is throwing an exception the best way to indicate command generation failure?
    command.toTry match {
      case Success(ic) => ic
      case Failure(e) => throw new Exception("Failed command instantiation", e)
    }
  }

  /**
    * Redirect the stdout and stderr to the appropriate files. While not necessary, mark the job as not receiving any
    * stdin by pointing it at /dev/null.
    *
    * If the `command` errors for some reason, put a "-1" into the rc file.
    */
  def redirectOutputs(command: String): String = {
    // > 128 is the cutoff for signal-induced process deaths such as might be observed with abort.
    // http://www.tldp.org/LDP/abs/html/exitcodes.html
    s"""$command < /dev/null || { rc=$$?; if [ "$$rc" -gt "128" ]; then echo $$rc; else echo -1; fi } > ${jobPaths.returnCode}"""
  }

  /** A tag that may be used for logging. */
  lazy val tag = s"${this.getClass.getSimpleName} [UUID(${workflowIdForLogging.shortString}):${jobDescriptor.key.tag}]"

  /**
    * When returns true, the `remoteStdErrPath` will be read. If contents of that path are non-empty, the job will fail.
    *
    * @return True if a non-empty `remoteStdErrPath` should fail the job.
    */
  lazy val failOnStdErr: Boolean = RuntimeAttributesValidation.extract(
    FailOnStderrValidation.instance, validatedRuntimeAttributes)

  /**
    * Returns the behavior for continuing on the return code, obtained by converting `returnCodeContents` to an Int.
    *
    * @return the behavior for continuing on the return code.
    */
  lazy val continueOnReturnCode: ContinueOnReturnCode = RuntimeAttributesValidation.extract(
    ContinueOnReturnCodeValidation.instance, validatedRuntimeAttributes)

  /**
    * Returns the max number of times that a failed job should be retried, obtained by converting `maxRetries` to an Int.
    */
  lazy val maxRetries: Int = RuntimeAttributesValidation.extract(
    MaxRetriesValidation.instance, validatedRuntimeAttributes)


  lazy val previousFailedRetries: Int = jobDescriptor.prefetchedKvStoreEntries.get(BackendLifecycleActorFactory.FailedRetryCountKey) match {
    case Some(KvPair(_,v)) => v.toInt
    case _ => 0
  }

  /**
    * Execute the job specified in the params. Should return a `StandardAsyncPendingExecutionHandle`, or a
    * `FailedExecutionHandle`.
    *
    * @return the execution handle for the job.
    */
  def execute(): ExecutionHandle = {
    throw new UnsupportedOperationException(s"Neither execute() nor executeAsync() implemented by $getClass")
  }

  /**
    * Async execute the job specified in the params. Should return a `StandardAsyncPendingExecutionHandle`, or a
    * `FailedExecutionHandle`.
    *
    * @return the execution handle for the job.
    */
  def executeAsync(): Future[ExecutionHandle] = Future.fromTry(Try(execute()))

  /**
    * Recovers the specified job id, or starts a new job. The default implementation simply calls execute().
    *
    * @param jobId The previously recorded job id.
    * @return the execution handle for the job.
    */
  def recover(jobId: StandardAsyncJob): ExecutionHandle = execute()

  /**
    * Async recovers the specified job id, or starts a new job. The default implementation simply calls execute().
    *
    * @param jobId The previously recorded job id.
    * @return the execution handle for the job.
    */
  def recoverAsync(jobId: StandardAsyncJob): Future[ExecutionHandle] = Future.fromTry(Try(recover(jobId)))

  /**
    * A correct implementation should reconnect to the specified job id, and upon reconnection try to abort the job.
    * The job might already be aborted (or in any other state) and this method should be robust to that scenario.
    * This method is needed in case Cromwell restarts while a workflow was being aborted. The engine cannot guarantee
    * that all backend actors will have time to receive and process the abort command before the server shuts down.
    * For that reason, upon restart, this method should always try to abort the job.
    *
    * The default implementation returns a JobReconnectionNotSupportedException failure which will result in a job failure.
    *
    * @param jobId The previously recorded job id.
    * @return the execution handle for the job.
    */
  def reconnectToAbortAsync(jobId: StandardAsyncJob): Future[ExecutionHandle] = Future.failed(JobReconnectionNotSupportedException(jobDescriptor.key))

  /**
    * This is in spirit similar to recover except it does not defaults back to running the job if not implemented.
    * This is used after a server restart to reconnect to jobs while the workflow was Failing (because another job failed). The workflow is bound
    * to fail eventually for that reason but in the meantime we want to reconnect to running jobs to update their status.
    *
    * The default implementation returns a JobReconnectionNotSupportedException failure which will result in a job failure.
    *
    * @param jobId The previously recorded job id.
    * @return the execution handle for the job.
    */
  def reconnectAsync(jobId: StandardAsyncJob): Future[ExecutionHandle] = Future.failed(JobReconnectionNotSupportedException(jobDescriptor.key))

  /**
    * Returns the run status for the job.
    *
    * @param handle The handle of the running job.
    * @return The status of the job.
    */
  def pollStatus(handle: StandardAsyncPendingExecutionHandle): StandardAsyncRunState = {
    throw new UnsupportedOperationException(s"Neither pollStatus nor pollStatusAsync implemented by $getClass")
  }

  /**
    * Returns the async run status for the job.
    *
    * @param handle The handle of the running job.
    * @return The status of the job.
    */
  def pollStatusAsync(handle: StandardAsyncPendingExecutionHandle): Future[StandardAsyncRunState] = Future.fromTry(Try(pollStatus(handle)))

  /**
    * Adds custom behavior invoked when polling fails due to some exception. By default adds nothing.
    *
    * Examples may be when certain error codes are detected during polling, and a specific handle should be returned.
    *
    * For example, say a custom JobNotFoundException should be mapped to a FailedNonRetryableExecutionHandle.
    *
    * @return A partial function handler for exceptions after polling.
    */
  def customPollStatusFailure: PartialFunction[(ExecutionHandle, Exception), ExecutionHandle] = {
    PartialFunction.empty
  }

  /**
    * Returns true when a job is complete, either successfully or unsuccessfully.
    *
    * @param runStatus The run status.
    * @return True if the job has completed.
    */
  def isTerminal(runStatus: StandardAsyncRunState): Boolean

  /**
    * Returns any events retrieved from the terminal run status.
    *
    * @param runStatus The terminal run status, as defined by isTerminal.
    * @return The execution events.
    */
  def getTerminalEvents(runStatus: StandardAsyncRunState): Seq[ExecutionEvent] = Seq.empty

  /**
    * Returns true if the status represents a completion.
    *
    * @param runStatus The run status.
    * @return True if the job is done.
    */
  def isDone(runStatus: StandardAsyncRunState): Boolean = true

  /**
    * Returns any custom metadata from the polled status.
    *
    * @param runStatus The run status.
    * @return The job metadata.
    */
  def getTerminalMetadata(runStatus: StandardAsyncRunState): Map[String, Any] = Map.empty

  /**
    * Attempts to abort a job when an abort signal is retrieved.
    *
    * If `requestsAbortAndDiesImmediately` is true, then the actor will die immediately after this method returns.
    *
    * @param jobId The job to abort.
    */
  def tryAbort(jobId: StandardAsyncJob): Unit = {}

  /**
    * Returns true if when an abort signal is retrieved, the actor makes an attempt to abort and then immediately stops
    * itself _without_ polling for an aborted status.
    *
    * NOTE: When this value is set to `false`, `tryAbort` must write the rc and stderr files, as they will still be
    * processed during a poll that returns a terminal status.
    *
    * The default is true.
    *
    * @return true if actor should request an abort and then die immediately.
    */
  def requestsAbortAndDiesImmediately: Boolean = true

  /**
    * Return true if the return code is an abort code.
    *
    * By default, return codes `SIGINT` and `SIGTERM` return true.
    *
    * @param returnCode The return code.
    * @return True if the return code is for an abort.
    */
  def isAbort(returnCode: Int): Boolean = returnCode == SIGINT || returnCode == SIGTERM || returnCode == SIGKILL

  /**
    * Custom behavior to run after an abort signal is processed.
    *
    * By default handles the behavior of `requestsAbortAndDiesImmediately`.
    */
  def postAbort(): Unit = {
    if (requestsAbortAndDiesImmediately) {
      tellMetadata(Map(CallMetadataKeys.BackendStatus -> "Aborted"))
      context.parent ! JobAbortedResponse(jobDescriptor.key)
      context.stop(self)
    }
  }

  /**
    * Output value mapper.
    *
    * @param womValue The original WOM value.
    * @return The Try wrapped and mapped WOM value.
    */
  final def outputValueMapper(womValue: WomValue): Try[WomValue] = {
    WomFileMapper.mapWomFiles(mapOutputWomFile, Set.empty)(womValue)
  }

  /**
    * Used to convert to output paths.
    *
    */
  def mapOutputWomFile(womFile: WomFile): WomFile =
    womFile.mapFile{
      path =>
        val pathFromContainerInputs = jobPaths.hostPathFromContainerInputs(path)
        pathFromContainerInputs.toAbsolutePath.toString
    }

  /**
    * Tries to evaluate the outputs.
    *
    * Used during handleExecutionSuccess.
    *
    * @return A Try wrapping evaluated outputs.
    */
  def evaluateOutputs()(implicit ec: ExecutionContext): Future[EvaluatedJobOutputs] = {
    OutputEvaluator.evaluateOutputs(jobDescriptor, backendEngineFunctions, outputValueMapper)
  }

  /**
    * Tests whether an attempted result of evaluateOutputs should possibly be retried.
    *
    * If the exception is a CromwellAggregatedException, this method will recursively call into itself checking if the
    * inner exceptions should be retried by using retryEvaluateOutputs.
    *
    * Custom implementations of handleExecutionSuccess should use this method when testing the results of
    * evaluateOutputs.
    *
    * However, to implement the actual check, override the function retryEvaluateOutputs.
    *
    * @param exception The exception, possibly an CromwellAggregatedException.
    * @return True if evaluateOutputs should be retried later.
    */
  final def retryEvaluateOutputsAggregated(exception: Exception): Boolean = {
    exception match {
      case aggregated: CromwellAggregatedException =>
        aggregated.throwables.collectFirst {
          case exception: Exception if retryEvaluateOutputsAggregated(exception) => exception
        }.isDefined
      case _ => retryEvaluateOutputs(exception)
    }
  }

  /**
    * Tests whether an attempted result of evaluateOutputs should possibly be retried.
    *
    * Override this function to check for different types of exceptions.
    *
    * Custom implementations of handleExecutionSuccess should NOT use method when testing the results of
    * evaluateOutputs, as this method does not recurse into aggregated exceptions. Instead custom overrides of
    * handleExecutionSuccess should call into retryEvaluateOutputsAggregated.
    *
    * @param exception The exception, possibly an internal instance retrieved from within a CromwellAggregatedException.
    * @return True if evaluateOutputs should be retried later.
    */
  def retryEvaluateOutputs(exception: Exception): Boolean = false

  /**
    * Process a successful run, as defined by `isSuccess`.
    *
    * @param runStatus  The run status.
    * @param handle     The execution handle.
    * @param returnCode The return code.
    * @return The execution handle.
    */
  def handleExecutionSuccess(runStatus: StandardAsyncRunState,
                             handle: StandardAsyncPendingExecutionHandle,
                             returnCode: Int)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    evaluateOutputs() map {
      case ValidJobOutputs(outputs) =>
        // Need to make sure the paths are up to date before sending the detritus back in the response
        updateJobPaths()
        SuccessfulExecutionHandle(outputs, returnCode, jobPaths.detritusPaths, getTerminalEvents(runStatus))
      case InvalidJobOutputs(errors) =>
        val exception = new MessageAggregation {
          override def exceptionContext: String = "Failed to evaluate job outputs"
          override def errorMessages: Traversable[String] = errors.toList
        }
        FailedNonRetryableExecutionHandle(exception)
      case JobOutputsEvaluationException(exception: Exception) if retryEvaluateOutputsAggregated(exception) =>
        // Return the execution handle in this case to retry the operation
        handle
      case JobOutputsEvaluationException(ex) => FailedNonRetryableExecutionHandle(ex)
    }
  }


  /**
    * Process an unsuccessful run, as interpreted by `handleExecutionFailure`.
    *
    * @param runStatus The run status.
    * @return The execution handle.
    */
  def retryElseFail(runStatus: StandardAsyncRunState,
                    backendExecutionStatus: Future[ExecutionHandle]): Future[ExecutionHandle] = {

    val retryable = previousFailedRetries < maxRetries

    backendExecutionStatus flatMap {
      case failed: FailedNonRetryableExecutionHandle if retryable =>
        incrementFailedRetryCount map { _ =>
          FailedRetryableExecutionHandle(failed.throwable, failed.returnCode)
        }
      case _ => backendExecutionStatus
    }
  }

  /**
    * Process an unsuccessful run, as defined by `isDone`.
    *
    * @param runStatus The run status.
    * @return The execution handle.
    */
  def handleExecutionFailure(runStatus: StandardAsyncRunState,
                             returnCode: Option[Int]): Future[ExecutionHandle] = {
    val exception = new RuntimeException(s"Task ${jobDescriptor.key.tag} failed for unknown reason: $runStatus")
    Future.successful(FailedNonRetryableExecutionHandle(exception, returnCode))
  }

  // See executeOrRecoverSuccess
  private var missedAbort = false
  private case class CheckMissedAbort(jobId: StandardAsyncJob)

  context.become(kvClientReceive orElse standardReceiveBehavior(None) orElse slowJobWarningReceive orElse receive)

  def standardReceiveBehavior(jobIdOption: Option[StandardAsyncJob]): Receive = LoggingReceive {
    case AbortJobCommand =>
      jobIdOption match {
        case Some(jobId) =>
          Try(tryAbort(jobId)) match {
            case Success(_) => jobLogger.info("{} Aborted {}", tag: Any, jobId)
            case Failure(ex) => jobLogger.warn("{} Failed to abort {}: {}", tag, jobId, ex.getMessage)
          }
        case None => missedAbort = true
      }
      postAbort()
    case CheckMissedAbort(jobId: StandardAsyncJob) =>
      context.become(kvClientReceive orElse standardReceiveBehavior(Option(jobId)) orElse receive)
      if (missedAbort)
        self ! AbortJobCommand
  }

  override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    val executeOrRecoverFuture = {
      mode match {
        case Reconnect(jobId: StandardAsyncJob@unchecked) => reconnectAsync(jobId)
        case ReconnectToAbort(jobId: StandardAsyncJob@unchecked) => reconnectToAbortAsync(jobId)
        case Recover(jobId: StandardAsyncJob@unchecked) => recoverAsync(jobId)
        case _ =>
          tellMetadata(startMetadataKeyValues)
          executeAsync()
      }
    }

    executeOrRecoverFuture flatMap executeOrRecoverSuccess recoverWith {
      case throwable: Throwable => Future failed {
        jobLogger.error(s"Error attempting to $mode", throwable)
        throwable
      }
    }
  }

  private def executeOrRecoverSuccess(executionHandle: ExecutionHandle): Future[ExecutionHandle] = {
    executionHandle match {
      case handle: PendingExecutionHandle[StandardAsyncJob@unchecked, StandardAsyncRunInfo@unchecked, StandardAsyncRunState@unchecked] =>

        configurationDescriptor.slowJobWarningAfter foreach { duration => self ! WarnAboutSlownessAfter(handle.pendingJob.jobId, duration) }

        tellKvJobId(handle.pendingJob) map { _ =>
          jobLogger.info(s"job id: ${handle.pendingJob.jobId}")
          tellMetadata(Map(CallMetadataKeys.JobId -> handle.pendingJob.jobId))
          /*
          NOTE: Because of the async nature of the Scala Futures, there is a point in time where we have submitted this or
          the prior runnable to the thread pool this actor doesn't know the job id for aborting. These runnables are
          queued up and may still be run by the thread pool anytime in the future. Issue #1218 may address this
          inconsistency at a later time. For now, just go back and check if we missed the abort command.
          */
          self ! CheckMissedAbort(handle.pendingJob)
          executionHandle
        }
      case _ => Future.successful(executionHandle)
    }
  }

  override def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    previous match {
      case handle: PendingExecutionHandle[
        StandardAsyncJob@unchecked, StandardAsyncRunInfo@unchecked, StandardAsyncRunState@unchecked] =>

        jobLogger.debug(s"$tag Polling Job ${handle.pendingJob}")
        pollStatusAsync(handle) flatMap {
          backendRunStatus =>
            self ! WarnAboutSlownessIfNecessary
            handlePollSuccess(handle, backendRunStatus)
        } recover {
          case throwable =>
            handlePollFailure(handle, throwable)
        }
      case successful: SuccessfulExecutionHandle => Future.successful(successful)
      case failed: FailedNonRetryableExecutionHandle => Future.successful(failed)
      case failedRetryable: FailedRetryableExecutionHandle => Future.successful(failedRetryable)
      case badHandle => Future.failed(new IllegalArgumentException(s"Unexpected execution handle: $badHandle"))
    }
  }

  /**
    * Process a poll success.
    *
    * @param oldHandle The previous execution handle.
    * @param state The updated run state.
    * @return The updated execution handle.
    */
  def handlePollSuccess(oldHandle: StandardAsyncPendingExecutionHandle,
                        state: StandardAsyncRunState): Future[ExecutionHandle] = {
    val previousState = oldHandle.previousState
    if (!(previousState exists statusEquivalentTo(state))) {
      // If this is the first time checking the status, we log the transition as '-' to 'currentStatus'. Otherwise just use
      // the state names.
      // This logging and metadata publishing assumes that StandardAsyncRunState subtypes `toString` nicely to state names.
      val prevStatusName = previousState.map(_.toString).getOrElse("-")
      jobLogger.info(s"Status change from $prevStatusName to $state")
      tellMetadata(Map(CallMetadataKeys.BackendStatus -> state))
    }

    state match {
      case _ if isTerminal(state) =>
        val metadata = getTerminalMetadata(state)
        tellMetadata(metadata)
        handleExecutionResult(state, oldHandle)
      case s => Future.successful(oldHandle.copy(previousState = Option(s))) // Copy the current handle with updated previous status.
    }
  }

  /**
    * Process a poll failure.
    *
    * @param oldHandle The previous execution handle.
    * @param throwable The cause of the polling failure.
    * @return The updated execution handle.
    */
  def handlePollFailure(oldHandle: StandardAsyncPendingExecutionHandle,
                        throwable: Throwable): ExecutionHandle = {
    throwable match {
      case exception: Exception =>
        val handler: PartialFunction[(ExecutionHandle, Exception), ExecutionHandle] =
          customPollStatusFailure orElse {
            case (_: ExecutionHandle, exception: Exception) if isFatal(exception) =>
              // Log exceptions and return the original handle to try again.
              jobLogger.warn(s"Fatal exception polling for status. Job will fail.", exception)
              FailedNonRetryableExecutionHandle(exception)
            case (handle: ExecutionHandle, exception: Exception) =>
              // Log exceptions and return the original handle to try again.
              jobLogger.warn(s"Caught non-fatal ${exception.getClass.getSimpleName} exception trying to poll, retrying", exception)
              handle
          }
        handler((oldHandle, exception))
      case error: Error => throw error // JVM-ending calamity.
      case _: Throwable =>
        // Someone has subclassed or instantiated Throwable directly. Kill the job. They should be using an Exception.
        FailedNonRetryableExecutionHandle(throwable)
    }
  }

  /**
    * Process an execution result.
    *
    * @param status The execution status.
    * @param oldHandle The previous execution handle.
    * @return The updated execution handle.
    */
  def handleExecutionResult(status: StandardAsyncRunState,
                            oldHandle: StandardAsyncPendingExecutionHandle): Future[ExecutionHandle] = {

    val stderr = jobPaths.standardPaths.error
    lazy val stderrAsOption: Option[Path] = Option(stderr)

    val stderrSizeAndReturnCode = for {
      returnCodeAsString <- asyncIo.contentAsStringAsync(jobPaths.returnCode, None, failOnOverflow = false)
      // Only check stderr size if we need to, otherwise this results in a lot of unnecessary I/O that
      // may fail due to race conditions on quickly-executing jobs.
      stderrSize <- if (failOnStdErr) asyncIo.sizeAsync(stderr) else Future.successful(0L)
    } yield (stderrSize, returnCodeAsString)

    stderrSizeAndReturnCode flatMap {
      case (stderrSize, returnCodeAsString) =>
        val tryReturnCodeAsInt = Try(returnCodeAsString.trim.toInt)

        if (isDone(status)) {
          tryReturnCodeAsInt match {
            case Success(returnCodeAsInt) if failOnStdErr && stderrSize.intValue > 0 =>
              val executionHandle = Future.successful(FailedNonRetryableExecutionHandle(StderrNonEmpty(jobDescriptor.key.tag, stderrSize, stderrAsOption), Option(returnCodeAsInt)))
              retryElseFail(status, executionHandle)
            case Success(returnCodeAsInt) if isAbort(returnCodeAsInt) =>
              Future.successful(AbortedExecutionHandle)
            case Success(returnCodeAsInt) if !continueOnReturnCode.continueFor(returnCodeAsInt) =>
              val executionHandle = Future.successful(FailedNonRetryableExecutionHandle(WrongReturnCode(jobDescriptor.key.tag, returnCodeAsInt, stderrAsOption), Option(returnCodeAsInt)))
              retryElseFail(status, executionHandle)
            case Success(returnCodeAsInt) =>
              handleExecutionSuccess(status, oldHandle, returnCodeAsInt)
            case Failure(_) =>
              Future.successful(FailedNonRetryableExecutionHandle(ReturnCodeIsNotAnInt(jobDescriptor.key.tag, returnCodeAsString, stderrAsOption)))
          }
        } else {
          val failureStatus = handleExecutionFailure(status, tryReturnCodeAsInt.toOption)
          retryElseFail(status, failureStatus)
        }
    } recoverWith {
      case exception =>
        if (isDone(status)) Future.successful(FailedNonRetryableExecutionHandle(exception))
        else {
          val failureStatus = handleExecutionFailure(status, None)
          retryElseFail(status, failureStatus)
        }
    }
  }

  /**
    * Send the job id of the running job to the key value store.
    *
    * @param runningJob The running job.
    */
  def tellKvJobId(runningJob: StandardAsyncJob): Future[KvResponse] = {
    val kvJobKey =
      KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt)
    val scopedKey = ScopedKey(jobDescriptor.workflowDescriptor.id, kvJobKey, jobIdKey)
    val kvValue = runningJob.jobId
    val kvPair = KvPair(scopedKey, kvValue)
    val kvPut = KvPut(kvPair)
    makeKvRequest(Seq(kvPut)).map(_.head)
  }

  /**
    * Increment the retry count for this failed job in the key value store.
    */
  def incrementFailedRetryCount: Future[KvResponse] = {
    val futureKvJobKey =
      KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt + 1)
    val scopedKey = ScopedKey(jobDescriptor.workflowDescriptor.id, futureKvJobKey, BackendLifecycleActorFactory.FailedRetryCountKey)
    val kvValue = (previousFailedRetries + 1).toString
    val kvPair = KvPair(scopedKey, kvValue)
    val kvPut = KvPut(kvPair)
    makeKvRequest(Seq(kvPut)).map(_.head)
  }

  /**
    * Sends metadata to the metadata store.
    *
    * @param metadataKeyValues Key/Values to store.
    */
  def tellMetadata(metadataKeyValues: Map[String, Any]): Unit = {
    import cromwell.services.metadata.MetadataService.implicits.MetadataAutoPutter
    serviceRegistryActor.putMetadata(jobDescriptor.workflowDescriptor.id, Option(jobDescriptor.key), metadataKeyValues)
  }

  override protected implicit lazy val ec: ExecutionContextExecutor = context.dispatcher
}

/**
  * Implements the marker trait for a job id using a String.
  *
  * @param jobId The job id.
  */
final case class StandardAsyncJob(jobId: String) extends JobId
