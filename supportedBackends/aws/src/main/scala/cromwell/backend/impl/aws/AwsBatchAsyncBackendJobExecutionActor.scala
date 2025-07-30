/*
 * Copyright 2018 Amazon.com, Inc. or its affiliates.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 *  BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 *  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 *  THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 *  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 *  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

package cromwell.backend.impl.aws

import java.net.{SocketTimeoutException, URLDecoder}
import java.io.FileNotFoundException
import java.nio.file.Paths
import akka.actor.ActorRef
import akka.pattern.AskSupport
import akka.util.Timeout
import cats.implicits._
import common.exception.MessageAggregation
import common.collections.EnhancedCollections._
import common.util.StringUtil._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cromwell.backend._
import cromwell.backend.async._
import cromwell.backend.impl.aws.IntervalLimitedAwsJobSubmitActor.SubmitAwsJobRequest
import cromwell.backend.impl.aws.OccasionalStatusPollingActor.{NotifyOfStatus, WhatsMyStatus}
import cromwell.backend.impl.aws.RunStatus.{Initializing, TerminalRunStatus}
import cromwell.backend.impl.aws.io._
import cromwell.backend.io.DirectoryFunctions
import cromwell.backend.io.JobPaths
import cromwell.backend.standard.{
  ScriptPreambleData,
  StandardAsyncExecutionActor,
  StandardAsyncExecutionActorParams,
  StandardAsyncJob
}
import cromwell.backend.OutputEvaluator._
import cromwell.backend.standard.retry.memory.MemoryRetryResult
import cromwell.core._
import cromwell.core.path.{DefaultPathBuilder, Path, PathBuilder, PathFactory}
import cromwell.core.io.{DefaultIoCommandBuilder, IoCommandBuilder}
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.filesystems.s3.S3Path
import cromwell.filesystems.s3.batch.S3BatchCommandBuilder
import cromwell.services.keyvalue.KvClient
import cromwell.services.metadata.CallMetadataKeys
import org.slf4j.{Logger, LoggerFactory}
import software.amazon.awssdk.services.batch.BatchClient
import software.amazon.awssdk.services.batch.model._
import wom.callable.Callable.OutputDefinition
import wom.callable.MetaValueElement.{MetaValueElementBoolean, MetaValueElementObject}
import wom.core.FullyQualifiedName
import wom.expression.NoIoFunctionSet
import wom.types.{
  WomArrayType,
  WomCompositeType,
  WomMapType,
  WomOptionalType,
  WomPairType,
  WomPrimitiveFileType,
  WomPrimitiveType,
  WomSingleFileType,
  WomType
}
import wom.values._

import scala.concurrent._
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.control.NoStackTrace
import scala.util.{Failure, Success, Try}

/**
  * The `AwsBatchAsyncBackendJobExecutionActor` creates and manages a job. The job itself is encapsulated by the
  * functionality in `AwsBatchJob`
  */
object AwsBatchAsyncBackendJobExecutionActor {
  val AwsBatchOperationIdKey = "__aws_batch_operation_id"

  type AwsBatchPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, AwsBatchJob, RunStatus]
}

class AwsBatchAsyncBackendJobExecutionActor(
  override val standardParams: StandardAsyncExecutionActorParams
) extends BackendJobLifecycleActor
    with StandardAsyncExecutionActor
    with AwsBatchJobCachingActorHelper
    with KvClient
    with AskSupport
    with AwsPlatform {

  /** The builder for `IoCommands` to the storage system used by jobs executed
   * by this backend
   */
  override lazy val ioCommandBuilder: IoCommandBuilder =
    configuration.fileSystem match {
      case AWSBatchStorageSystems.s3 => S3BatchCommandBuilder
      case _ => DefaultIoCommandBuilder
    }

  // the cromwell backend Actor
  val backendSingletonActor: ActorRef =
    standardParams.backendSingletonActorOption.getOrElse(
      throw new RuntimeException(
        s"AWS Backend actor cannot exist without its backend singleton (of type ${AwsBatchSingletonActor.getClass.getSimpleName})"
      )
    )

  import AwsBatchAsyncBackendJobExecutionActor._

  val Log: Logger = LoggerFactory.getLogger(AwsBatchAsyncBackendJobExecutionActor.getClass)

  override type StandardAsyncRunInfo = AwsBatchJob

  override type StandardAsyncRunState = RunStatus

  /**
   * Determines if two run statuses are equal
   *
   * @param thiz a `RunStatus`
   * @param that a `RunStatus`
   * @return true if they are equal, else false
   */
  def statusEquivalentTo(thiz: StandardAsyncRunState)(that: StandardAsyncRunState): Boolean = thiz == that

  override lazy val pollBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(1.second, 5.minutes, 1.1)

  override lazy val executeOrRecoverBackOff: SimpleExponentialBackoff =
    SimpleExponentialBackoff(initialInterval = 3 seconds, maxInterval = 20 seconds, multiplier = 1.1)

  // the name (String) of the docker image that will be used to contain this job
  private lazy val jobDockerImage =
    jobDescriptor.maybeCallCachingEligible.dockerHash.getOrElse(runtimeAttributes.dockerImage)

  override lazy val dockerImageUsed: Option[String] = Option(jobDockerImage)

  // |cd ${jobPaths.script.parent.pathWithoutScheme}; ls | grep -v script | xargs rm -rf; cd -

  // overriden function : used in StandardAsyncExecutionActor
  override def inputsToNotLocalize: Set[WomFile] =
    jobDescriptor.fullyQualifiedInputs.collect {
      case (_, womFile: WomFile)
          if jobDescriptor
            .findInputFilesByParameterMeta {
              case MetaValueElementObject(values) =>
                values.get("localization_optional").contains(MetaValueElementBoolean(true))
              case _ => false
            }
            .contains(womFile) =>
        womFile

    }.toSet

  // custom routine that simply returns a boolean. Used here in callInputFiles
  private def isLocalizationOptional: WomValue => Boolean =
    jobDescriptor.fullyQualifiedInputs.collect {
      case (_, womFile: WomFile)
          if jobDescriptor
            .findInputFilesByParameterMeta {
              case MetaValueElementObject(values) =>
                values.get("localization_optional").contains(MetaValueElementBoolean(true))
              case _ => false
            }
            .contains(womFile) =>
        womFile

    }.toSet

  private lazy val execScript =
    s"""|#!$jobShell
        |find ${jobPaths.script.parent.pathWithoutScheme} -group root | grep -v script | xargs rm -vrf
        |${jobPaths.script.pathWithoutScheme}
        |""".stripMargin

  /* Batch job object (see AwsBatchJob). This has the configuration necessary
   * to perform all operations with the AWS Batch infrastructure. This is
   * where the real work happens
   *
   * Rundown of the command string:
   *
   * commandScriptContents: This is ErrorOr[String] that includes a full
   *                        bash script designed to get output into a file
   *                        It's defined in JobPaths.scala. Other backends
   *                        do this funky "write to a file in the storage service,
   *                        have the container pick up that file and run it" thing.
   *
   *                        But, I'm not convinced yet that Cromwell needs this,
   *                        and I think that we can pass over to AWS Batch
   *                        what is needed to run. So...why do anything like
   *                        the following?
   *
   * commandScriptContents.fold(
   *   errors => Future.failed(new RuntimeException(errors.toList.mkString(", "))),
   *   callPaths.script.write)
   *
   *
   * So, what we're passing to AwsBatchJob here is the literal command string -
   *
   * instantiatedCommand.commandString: This is an InstantiatedCommand class and holds
   *                                    all the things about the command. It's defined
   *                                    in StandardAsyncExecutionActor
   *
   * NOTE: In order to get output from the command, the commandScriptContents
   * needs to push stuff out to S3. This is why we will eventually need
   * commandScriptContents here
   */

  /* part of the full commandScriptContents is overriden here, in the context of mixed S3/EFS support with globbing.
      we'll see how much we need...
   */

  lazy val cmdScript = configuration.fileSystem match {
    case AWSBatchStorageSystems.s3 => commandScriptContents.toEither.toOption.get
    case _ => execScript
  }

  lazy val batchJob: AwsBatchJob =
    AwsBatchJob(
      jobDescriptor,
      runtimeAttributes,
      instantiatedCommand.commandString,
      cmdScript,
      rcPath.toString,
      executionStdout,
      executionStderr,
      generateAwsBatchInputs(jobDescriptor),
      generateAwsBatchOutputs(jobDescriptor),
      jobPaths,
      Seq.empty[AwsBatchParameter],
      configuration.awsConfig.region,
      Option(configuration.awsAuth),
      configuration.fsxMntPoint,
      configuration.efsMntPoint,
      Option(runtimeAttributes.efsMakeMD5),
      Option(runtimeAttributes.efsDelocalize),
      Option(runtimeAttributes.tagResources),
      runtimeAttributes.logGroupName,
      runtimeAttributes.additionalTags
    )

  // setup batch client to query job container info
  lazy val batchClient: BatchClient = {
    val builder = BatchClient.builder()
    configureClient(builder, batchJob.optAwsAuthMode, batchJob.configRegion)
  }

  /* Tries to abort the job in flight
   *
   * @param job A StandardAsyncJob object (has jobId value) to cancel
   * @return Nothing
   *
   */
  override def tryAbort(job: StandardAsyncJob): Unit = {
    batchJob.abort(
      job.jobId
    ) // job.JobId should be the AWS Batch Job Id based on analysis of other backends
    Log.info(
      s"Attempted CancelJob operation in AWS Batch for Job ID ${job.jobId}. There were no errors during the operation"
    )
    Log.info(
      s"We have normality. Anything you still can't cope with is therefore your own problem"
    )
    Log.info(s"https://www.youtube.com/watch?v=YCRxnjE7JVs")
    ()
  }

  override def requestsAbortAndDiesImmediately: Boolean = false

  /**
   * Takes two arrays of remote and local WOM File paths and generates the necessary AwsBatchInputs.
   */

  // taken from the WomExpression trait:
  private def areAllFilesOptional(womType: WomType): Boolean = {
    def innerAreAllFileTypesInWomTypeOptional(womType: WomType): Boolean = womType match {
      case WomOptionalType(_: WomPrimitiveFileType) =>
        true
      case _: WomPrimitiveFileType =>
        false
      case _: WomPrimitiveType =>
        true // WomPairTypes and WomCompositeTypes may have non-File components here which is fine.
      case WomArrayType(inner) => innerAreAllFileTypesInWomTypeOptional(inner)
      case WomMapType(_, inner) => innerAreAllFileTypesInWomTypeOptional(inner)
      case WomPairType(leftType, rightType) =>
        innerAreAllFileTypesInWomTypeOptional(leftType) && innerAreAllFileTypesInWomTypeOptional(rightType)
      case WomCompositeType(typeMap, _) => typeMap.values.forall(innerAreAllFileTypesInWomTypeOptional)
      case _ => false
    }

    // At the outermost level, primitives are never optional.
    womType match {
      case _: WomPrimitiveType =>
        false
      case _ => innerAreAllFileTypesInWomTypeOptional(womType)
    }
  }

  private def inputsFromWomFiles(
    namePrefix: String,
    remotePathArray: Seq[WomFile],
    localPathArray: Seq[WomFile],
    isOptional: Seq[Boolean], // can file be missing (File?)
    isLocOptional: Seq[Boolean] // can localization be skipped (through meta tags)
  ): Iterable[AwsBatchInput] =
    (remotePathArray zip localPathArray zip isOptional zip isLocOptional zipWithIndex) flatMap {
      case ((((remotePath, localPath), optional), locOptional), index) =>
        var localPathString = localPath.valueString
        // not localizing : keep remote path
        if (locOptional) {
          localPathString = remotePath.valueString
          // localizing : strip s3 prefixes
        } else if (localPathString.startsWith("s3://")) {
          localPathString = localPathString.replace("s3://", "")
        } else if (localPathString.startsWith("s3:/")) {
          localPathString = localPathString.replace("s3:/", "")
        }
        val localBuiltPath = DefaultPathBuilder.get(localPathString)
        Seq(
          AwsBatchFileInput(
            s"$namePrefix-$index",
            remotePath.valueString,
            localBuiltPath,
            workingDisk,
            optional,
            locOptional
          )
        )
    }

  /**
    * Turns WomFiles into relative paths.  These paths are relative to the working disk.
    *
    * relativeLocalizationPath("foo/bar.txt") -> "foo/bar.txt"
    * relativeLocalizationPath("s3://some/bucket/foo.txt") -> "some/bucket/foo.txt"
    */
  override protected def relativeLocalizationPath(file: WomFile): WomFile =
    file.mapFile(value =>
      getPath(value) match {
        // for s3 paths :
        case Success(path: S3Path) =>
          configuration.fileSystem match {
            case AWSBatchStorageSystems.s3 =>
              URLDecoder.decode(path.pathWithoutScheme, "UTF-8")
            case _ =>
              URLDecoder.decode(path.toString, "UTF-8")
          }
        // non-s3 paths
        case _ =>
          URLDecoder.decode(value, "UTF-8")
      }
    )

  /**
    * Generate a set of inputs based on a job description
    * @param jobDescriptor the job descriptor from Cromwell
    * @return the inputs derived from the descriptor
    */
  private[aws] def generateAwsBatchInputs(
    jobDescriptor: BackendJobDescriptor
  ): Set[AwsBatchInput] = {
    val writeFunctionFiles = instantiatedCommand.createdFiles map { f =>
      f.file.value.md5SumShort -> List(f)
    } toMap

    val writeFunctionInputs = writeFunctionFiles flatMap { case (name, files) =>
      inputsFromWomFiles(
        name,
        files.map(_.file),
        files.map(localizationPath),
        Seq.fill(files.size)(false), // Add Seq() of false for each file (isOptional)
        Seq.fill(files.size)(false) // Add Seq() of false for each file (locOptional)
      )
    }

    val callInputFiles: Map[FullyQualifiedName, Seq[(WomFile, Boolean, Boolean)]] =
      jobDescriptor.fullyQualifiedInputs safeMapValues {
        // consider womfiles, optional womfiles, lists of womfiles, optional lists of womfiles, lists of optional womfiles
        case womFile =>
          // we can skip optional_localization files here, but that prevents checking in the call script for existence.
          //   => pass on all input files, check optional_localization during the actual localization.
          val arrays: Seq[WomArray] = womFile collectAsSeq {
            case womFile: WomFile =>
              val files: List[WomSingleFile] = DirectoryFunctions
                .listWomSingleFiles(womFile, callPaths.workflowPaths)
                .toTry(s"Error getting single files for $womFile")
                .get
              WomArray(WomArrayType(WomSingleFileType), files)
            case _ =>
              WomArray(WomArrayType(WomSingleFileType), List.empty)
          }
          // if files found: define optional statuses
          val (isOptional, isLocOptional) = arrays.nonEmpty match {
            case true =>
              // mixed optional is cast to mandatory
              val isOptional = areAllFilesOptional(womFile.womType)
              val isLocOptional = isLocalizationOptional(womFile)
              (isOptional, isLocOptional)
            case false =>
              (false, false)
          }

          // Flatten and collect files along with the outer womFile's optional status and localization requirement
          arrays.flatMap(_.value).collect { case file: WomFile =>
            (file, isOptional, isLocOptional)
          }

      }
    val callInputInputs = callInputFiles flatMap { case (name, filesWithOptional) =>
      val files = filesWithOptional.map(_._1) // Extract the WomFiles
      val isOptional = filesWithOptional.map(_._2) // Extract the corresponding optional status
      val isLocOptional = filesWithOptional.map(_._3) // Extract the corresponding localization status

      inputsFromWomFiles(
        name,
        files,
        files.map(relativeLocalizationPath),
        isOptional, // Pass optional status
        isLocOptional // Pass localization requirement
      )
    }
    // this is a list : AwsBatchInput(name_in_wf, origin_such_as_s3, target_in_docker_relative, target_in_docker_disk[name mount], is_optional_file )
    val scriptInput: AwsBatchInput = AwsBatchFileInput(
      "script",
      jobPaths.script.pathAsString,
      DefaultPathBuilder.get(jobPaths.script.pathWithoutScheme),
      workingDisk,
      false,
      false
    )

    Set(scriptInput) ++ writeFunctionInputs ++ callInputInputs
  }

  /**
    * Given a path (relative or absolute), returns a (Path, AwsBatchVolume) tuple where the Path is
    * relative to the Volume's mount point
    *
    * @throws Exception if the `path` does not live in one of the supplied `disks`
    */
  private def relativePathAndVolume(
    path: String,
    disks: Seq[AwsBatchVolume]
  ): (Path, AwsBatchVolume) = {

    def getAbsolutePath(path: Path) =
      configuration.fileSystem match {
        case AWSBatchStorageSystems.s3 => AwsBatchWorkingDisk.MountPoint.resolve(path)
        case _ => AwsBatchWorkingDisk.MountPoint.resolve(path)
      }

    val absolutePath = DefaultPathBuilder.get(path) match {
      case p if !p.isAbsolute => getAbsolutePath(p)
      case p => p
    }
    disks.find(d => absolutePath.startsWith(d.mountPoint)) match {
      case Some(disk) => (disk.mountPoint.relativize(absolutePath), disk)
      case None =>
        throw new Exception(
          s"Absolute path $path doesn't appear to be under any mount points: ${disks.map(_.toString).mkString(", ")}"
        )
    }
  }

  /**
    * Produces names with a length less than 128 characters possibly by producing a digest of the name
    * @param referenceName the name to make safe
    * @return the name or the MD5sum of that name if the name is >= 128 characters
    */
  private def makeSafeAwsBatchReferenceName(referenceName: String) =
    if (referenceName.length <= 127) referenceName else referenceName.md5Sum

  private[aws] def generateAwsBatchOutputs(
    jobDescriptor: BackendJobDescriptor
  ): Set[AwsBatchFileOutput] = {
    import cats.syntax.validated._

    def evaluateFiles(output: OutputDefinition): List[(WomFile, Boolean)] = {
      // mixed mandatory/optional types are cast to mandatory.
      val is_optional = areAllFilesOptional(output.womType)

      Try(
        output.expression
          .evaluateFiles(
            jobDescriptor.localInputs,
            NoIoFunctionSet,
            output.womType
          )
          .map(_.toList map { womFile =>
            // Pair each WomFile with the optional status
            (womFile.file, is_optional)
          })
      ).getOrElse(List.empty[(WomFile, Boolean)].validNel)
        .getOrElse(List.empty)
    }

    val womFileOutputsEvaluated = jobDescriptor.taskCall.callable.outputs
      .flatMap(evaluateFiles)

    val womFileOutputs = womFileOutputsEvaluated map { case (file, isOptional) =>
      (relativeLocalizationPath(file), isOptional)
    }

    val flatmapped = womFileOutputs.distinct flatMap { case (file, isOptional) =>
      file.flattenFiles.map((_, isOptional))
    }
    // Generate AwsBatchFileOutput for each file based on type
    val outputs: Seq[AwsBatchFileOutput] = flatmapped flatMap {

      case (unlistedDirectory: WomUnlistedDirectory, _) =>
        generateUnlistedDirectoryOutputs(unlistedDirectory)
      case (singleFile: WomSingleFile, isOptional) =>
        generateAwsBatchSingleFileOutputs(singleFile, isOptional)
      case (globFile: WomGlobFile, _) =>
        generateAwsBatchGlobFileOutputs(globFile)
    }

    val additionalGlobOutput =
      jobDescriptor.taskCall.callable.additionalGlob.toList
        .flatMap(generateAwsBatchGlobFileOutputs)
        .toSet

    outputs.toSet ++ additionalGlobOutput
  }

  // used by generateAwsBatchOutputs, could potentially move this def within that function
  private def generateUnlistedDirectoryOutputs(womFile: WomUnlistedDirectory): List[AwsBatchFileOutput] = {
    val directoryPath = womFile.value.ensureSlashed
    val directoryListFile = womFile.value.ensureUnslashed + ".list"
    val dirDestinationPath = callRootPath.resolve(directoryPath).pathAsString
    val listDestinationPath = callRootPath.resolve(directoryListFile).pathAsString
    val (_, directoryDisk) = relativePathAndVolume(womFile.value, runtimeAttributes.disks)

    // We need both the collection directory and the collection list:
    List(
      // The collection directory:
      AwsBatchFileOutput(
        makeSafeAwsBatchReferenceName(directoryListFile),
        listDestinationPath,
        DefaultPathBuilder.get(directoryListFile),
        directoryDisk,
        false
      ),
      // The collection list file:
      AwsBatchFileOutput(
        makeSafeAwsBatchReferenceName(directoryPath),
        dirDestinationPath,
        DefaultPathBuilder.get(directoryPath + "*"),
        directoryDisk,
        false
      )
    )
  }

  // used by generateAwsBatchOutputs, could potentially move this def within that function
  private def generateAwsBatchSingleFileOutputs(
    womFile: WomSingleFile,
    isOptional: Boolean
  ): List[AwsBatchFileOutput] = {
    // rewrite this to create more flexibility
    //
    val destination = configuration.fileSystem match {
      case AWSBatchStorageSystems.s3 =>
        callRootPath.resolve(womFile.value.stripPrefix("/")).pathAsString
      case _ =>
        DefaultPathBuilder.get(womFile.valueString) match {
          case p if !p.isAbsolute =>
            callRootPath.resolve(womFile.value.stripPrefix("/")).pathAsString
          case p => p.pathAsString
        }

    }
    val (relpath, disk) = relativePathAndVolume(womFile.value, runtimeAttributes.disks)

    val output =
      if (
        configuration.efsMntPoint.isDefined &&
        configuration.efsMntPoint.getOrElse("").equals(disk.toString.split(" ")(1)) &&
        !runtimeAttributes.efsDelocalize
      ) {
        // name: String, s3key: String, local: Path, mount: AwsBatchVolume, optionalFile
        AwsBatchFileOutput(makeSafeAwsBatchReferenceName(womFile.value), womFile.value, relpath, disk, isOptional)
      } else {
        // if efs is not enabled, OR efs delocalization IS enabled, keep the s3 path as destination.
        AwsBatchFileOutput(makeSafeAwsBatchReferenceName(womFile.value),
                           URLDecoder.decode(destination, "UTF-8"),
                           relpath,
                           disk,
                           isOptional
        )
      }
    List(output)
  }

  // get a unique glob name locations & paths.
  //  1. globName :md5 hash of local PATH and WF_ID
  //  2. globbedDir : local path of the directory being globbed.
  //  3. volume on which the globbed data is located (eg root, efs, ...)
  //  4. target path for delocalization for globDir
  //  5. target path for delocalization for globList
  private def generateGlobPaths(womFile: WomGlobFile): (String, String, AwsBatchVolume, String, String) = {
    // add workflow id to hash for better conflict prevention
    val wfid = standardParams.jobDescriptor.toString.split(":")(0)
    val globName = GlobFunctions.globName(s"${womFile.value}-${wfid}")

    var globbedDirPath = Option(Paths.get(womFile.value).getParent()) match {
      case None => Paths.get(".")
      case Some(parent) => parent
    }
    while (globbedDirPath.toString().contains("*"))
      globbedDirPath = Option(globbedDirPath.getParent()) match {
        case None => Paths.get(".")
        case Some(parent) => parent
      }
    val globbedDir: String = globbedDirPath.toString()
    // generalize folder and list file
    val globDirectory = DefaultPathBuilder.get(globbedDir + "/." + globName + "/")
    val globListFile = DefaultPathBuilder.get(globbedDir + "/." + globName + ".list")

    // locate the disk where the globbed data resides
    val (_, globDirectoryDisk) = relativePathAndVolume(womFile.value, runtimeAttributes.disks)

    val (globDirectoryDestinationPath, globListFileDestinationPath) =
      if (
        configuration.efsMntPoint.isDefined &&
        configuration.efsMntPoint
          .getOrElse("")
          .equals(globDirectoryDisk.toString.split(" ")(1)) &&
        !runtimeAttributes.efsDelocalize
      ) {
        (globDirectory, globListFile)
      } else {
        // cannot resolve absolute paths : strip the leading '/'
        (
          callRootPath
            // first strip './' (glob in working dir), then '/' (relative path)
            .resolve(globDirectory.toString.stripPrefix("./").stripPrefix("/"))
            .pathAsString,
          callRootPath
            .resolve(globListFile.toString.stripPrefix("./").stripPrefix("/"))
            .pathAsString
        )
      }
    // return results
    (
      globName,
      globbedDir,
      globDirectoryDisk,
      globDirectoryDestinationPath.toString,
      globListFileDestinationPath.toString
    )

  }
  // used by generateAwsBatchOutputs, could potentially move this def within that function
  private def generateAwsBatchGlobFileOutputs(womFile: WomGlobFile): List[AwsBatchFileOutput] = {

    val (globName, globbedDir, globDirectoryDisk, globDirectoryDestinationPath, globListFileDestinationPath) =
      generateGlobPaths(womFile)
    val (relpathDir, _) = relativePathAndVolume(
      DefaultPathBuilder.get(globbedDir + "/." + globName + "/" + "*").toString,
      runtimeAttributes.disks
    )
    val (relpathList, _) = relativePathAndVolume(
      DefaultPathBuilder.get(globbedDir + "/." + globName + ".list").toString,
      runtimeAttributes.disks
    )
    // We need both the glob directory and the glob list:
    List(
      // The glob directory:
      AwsBatchFileOutput(DefaultPathBuilder.get(globbedDir.toString + "/." + globName + "/" + "*").toString,
                         globDirectoryDestinationPath,
                         relpathDir,
                         globDirectoryDisk,
                         false
      ),
      // The glob list file:
      AwsBatchFileOutput(DefaultPathBuilder.get(globbedDir.toString + "/." + globName + ".list").toString,
                         globListFileDestinationPath,
                         relpathList,
                         globDirectoryDisk,
                         false
      )
      // TODO: EVALUATE above vs below  (mainly the makeSafeAwsBatchReferenceName() routine)
      // The glob directory:
      // AwsBatchFileOutput(makeSafeAwsBatchReferenceName(globDirectory),
      //                    globDirectoryDestinationPath,
      //                    DefaultPathBuilder.get(globDirectory + "*"),
      //                    globDirectoryDisk
      // ),
      // The glob list file:
      // AwsBatchFileOutput(makeSafeAwsBatchReferenceName(globListFile),
      //                    globListFileDestinationPath,
      //                    DefaultPathBuilder.get(globListFile),
      //                    globDirectoryDisk
      // )
    )
  }

  override lazy val commandDirectory: Path = configuration.fileSystem match {
    case AWSBatchStorageSystems.s3 => AwsBatchWorkingDisk.MountPoint
    case _ => jobPaths.callExecutionRoot
  }

  override def scriptPreamble: ErrorOr[ScriptPreambleData] =
    configuration.fileSystem match {
      case AWSBatchStorageSystems.s3 => ScriptPreambleData("").validNel
      case _ => ScriptPreambleData("").validNel
    }

  override def scriptClosure: String =
    configuration.fileSystem match {
      case AWSBatchStorageSystems.s3 => ""
      case _ => s"exit $$(head -n 1 $rcPath)"
    }

  override def globParentDirectory(womGlobFile: WomGlobFile): Path =
    configuration.fileSystem match {
      case AWSBatchStorageSystems.s3 =>
        val (_, disk) = relativePathAndVolume(womGlobFile.value, runtimeAttributes.disks)
        disk.mountPoint
      case _ => commandDirectory
    }

  override def isTerminal(runStatus: RunStatus): Boolean =
    runStatus match {
      case _: TerminalRunStatus => true
      case _ => false
    }

  /**
    * Asynchronously upload the command script to the script path
    * @return a `Future` for the asynch operation
    */
  def uploadScriptFile(): Future[Unit] =
    commandScriptContents.fold(
      errors => Future.failed(new RuntimeException(errors.toList.mkString(", "))),
      asyncIo.writeAsync(jobPaths.script, _, Seq.empty)
    )

  // Primary entry point for cromwell to actually run something
  override def executeAsync(): Future[ExecutionHandle] =
    for {
      // upload the command script
      _ <- uploadScriptFile()
      completionPromise = Promise[SubmitJobResponse]()
      // send a message to the Actor requesting a job submission
      _ = backendSingletonActor ! SubmitAwsJobRequest(batchJob, attributes, completionPromise)
      // the future response of the submit job request
      submitJobResponse <- completionPromise.future
      // send a notify of status method to the Actor
      _ = backendSingletonActor ! NotifyOfStatus(runtimeAttributes.queueArn,
                                                 submitJobResponse.jobId,
                                                 Option(Initializing)
      )
    } yield PendingExecutionHandle(jobDescriptor,
                                   StandardAsyncJob(submitJobResponse.jobId),
                                   Option(batchJob),
                                   previousState = None
    )

  override def recoverAsync(jobId: StandardAsyncJob): Future[ExecutionHandle] = reconnectAsync(jobId)

  override def reconnectAsync(
    jobId: StandardAsyncJob
  ): Future[ExecutionHandle] = {
    val handle = PendingExecutionHandle[
      StandardAsyncJob,
      StandardAsyncRunInfo,
      StandardAsyncRunState
    ](jobDescriptor, jobId, Option(batchJob), previousState = None)
    Future.successful(handle)
  }

  override def reconnectToAbortAsync(jobId: StandardAsyncJob): Future[ExecutionHandle] = {
    tryAbort(jobId)
    reconnectAsync(jobId)
  }

  // This is called by Cromwell after initial execution (see executeAsync above)
  // It expects a Future[RunStatus]. In this case we'll simply call the
  // AWS Batch API to do this. The AwsBatchJob object in the PendingExecutionHandle
  // will have the actual status call, so our job here is simply to pull the
  // object out of the handle and execute the underlying status method
  override def pollStatusAsync(handle: AwsBatchPendingExecutionHandle): Future[RunStatus] = {
    val jobId = handle.pendingJob.jobId
    val job = handle.runInfo match {
      case Some(actualJob) => actualJob
      case None =>
        throw new RuntimeException(
          s"pollStatusAsync called but job not available. This should not happen. Job Id $jobId"
        )
    }

    implicit val timeout: Timeout = Timeout(5.seconds)

    def useQuickAnswerOrFallback(quick: Any): Future[RunStatus] = quick match {
      case NotifyOfStatus(_, _, Some(value)) =>
        Future.successful(value)
      case NotifyOfStatus(_, _, None) =>
        jobLogger.debug("Having to fall back to AWS query for status")
        Future.fromTry(job.status(jobId))
      case other =>
        val message =
          s"Programmer Error (please report this): Received an unexpected message from the OccasionalPollingActor: $other"
        jobLogger.error(message)
        Future.failed(new Exception(message) with NoStackTrace)
    }

    for {
      quickAnswer <- backendSingletonActor ? WhatsMyStatus(runtimeAttributes.queueArn, jobId)
      guaranteedAnswer <- useQuickAnswerOrFallback(quickAnswer)
    } yield guaranteedAnswer
  }

  override def handleExecutionResult(
    status: StandardAsyncRunState,
    oldHandle: StandardAsyncPendingExecutionHandle
  ): Future[ExecutionHandle] = {

    // get path to sderr
    val stderr = jobPaths.standardPaths.error
    lazy val stderrAsOption: Option[Path] = Option(stderr)
    // get the three needed variables, using helper functions below, or direct assignment.
    val stderrSizeAndReturnCodeAndMemoryRetry = for {
      returnCodeAsString <- JobExitCode
      // Only check stderr size if we need to, otherwise this results in a lot of unnecessary I/O that
      // may fail due to race conditions on quickly-executing jobs.
      stderrSize <- if (failOnStdErr) asyncIo.sizeAsync(stderr) else Future.successful(0L)
      retryWithMoreMemory <- memoryRetryRC(oldHandle.pendingJob)
    } yield (stderrSize, returnCodeAsString, retryWithMoreMemory)

    stderrSizeAndReturnCodeAndMemoryRetry flatMap { case (stderrSize, returnCodeAsString, retryWithMoreMemory) =>
      val tryReturnCodeAsInt = Try(returnCodeAsString.trim.toInt)
      jobLogger.debug(
        s"Handling execution Result with status '${status.toString()}' and returnCode ${returnCodeAsString}"
      )
      if (isDone(status)) {
        tryReturnCodeAsInt match {
          // stderr not empty : retry
          case Success(returnCodeAsInt) if failOnStdErr && stderrSize.intValue > 0 =>
            val executionHandle = Future.successful(
              FailedNonRetryableExecutionHandle(
                StderrNonEmpty(
                  jobDescriptor.key.tag,
                  stderrSize,
                  stderrAsOption
                ),
                Option(returnCodeAsInt),
                None
              )
            )
            retryElseFail(executionHandle)
          // job was aborted (cancelled by user?)
          // on AWS OOM kill are code 137 : check retryWithMoreMemory here
          case Success(returnCodeAsInt) if isAbort(returnCodeAsInt) && !retryWithMoreMemory =>
            jobLogger.debug(
              s"Job was aborted, code was : '${returnCodeAsString.stripLineEnd}'"
            )
            Future.successful(AbortedExecutionHandle)
          // if instance killed after RC.txt creation : edge case with status == Failed AND returnCode == [accepted values] => retry.
          case Success(returnCodeAsInt)
              if status.toString() == "Failed" && continueOnReturnCode
                .continueFor(returnCodeAsInt) =>
            jobLogger.debug(s"Suspected spot kill due to status/RC mismatch")
            val executionHandle = Future.successful(
              FailedNonRetryableExecutionHandle(
                UnExpectedStatus(
                  jobDescriptor.key.tag,
                  returnCodeAsInt,
                  status.toString(),
                  stderrAsOption
                ),
                Option(returnCodeAsInt),
                None
              )
            )
            retryElseFail(executionHandle)
          // job considered ok by accepted exit code
          case Success(returnCodeAsInt) if continueOnReturnCode.continueFor(returnCodeAsInt) =>
            handleExecutionSuccess(status, oldHandle, returnCodeAsInt)
          // job failed on out-of-memory : retry
          case Success(returnCodeAsInt) if retryWithMoreMemory =>
            jobLogger.info(
              s"Retrying job due to OOM with exit code : '${returnCodeAsString.stripLineEnd}' "
            )
            val executionHandle = Future.successful(
              FailedNonRetryableExecutionHandle(
                RetryWithMoreMemory(
                  jobDescriptor.key.tag,
                  stderrAsOption,
                  memoryRetryErrorKeys,
                  log
                ),
                Option(returnCodeAsInt),
                None
              )
            )
            retryElseFail(executionHandle,
                          MemoryRetryResult(retryWithMoreMemory, memoryRetryFactor, previousMemoryMultiplier)
            )
          // unaccepted return code : retry.
          case Success(returnCodeAsInt) =>
            jobLogger.debug(
              s"Retrying with wrong exit code : '${returnCodeAsString.stripLineEnd}'"
            )
            val executionHandle = Future.successful(
              FailedNonRetryableExecutionHandle(
                WrongReturnCode(
                  jobDescriptor.key.tag,
                  returnCodeAsInt,
                  stderrAsOption
                ),
                Option(returnCodeAsInt),
                None
              )
            )
            retryElseFail(executionHandle)
          case Failure(_) =>
            jobLogger.warn(
              s"General failure of job with exit code : '${returnCodeAsString.stripLineEnd}'"
            )
            Future.successful(
              FailedNonRetryableExecutionHandle(
                ReturnCodeIsNotAnInt(
                  jobDescriptor.key.tag,
                  returnCodeAsString,
                  stderrAsOption
                ),
                kvPairsToSave = None
              )
            )
        }
      } else {
        tryReturnCodeAsInt match {
          case Success(returnCodeAsInt) if retryWithMoreMemory && !continueOnReturnCode.continueFor(returnCodeAsInt) =>
            jobLogger.debug(s"job not done but retrying already? : ${status.toString()}")
            val executionHandle = Future.successful(
              FailedNonRetryableExecutionHandle(
                RetryWithMoreMemory(jobDescriptor.key.tag, stderrAsOption, memoryRetryErrorKeys, log),
                Option(returnCodeAsInt),
                None
              )
            )
            retryElseFail(executionHandle,
                          MemoryRetryResult(retryWithMoreMemory, memoryRetryFactor, previousMemoryMultiplier)
            )
          case _ =>
            val failureStatus = handleExecutionFailure(status, tryReturnCodeAsInt.toOption)
            retryElseFail(failureStatus)
        }
      }
    } recoverWith { case exception =>
      if (isDone(status)) Future.successful(FailedNonRetryableExecutionHandle(exception, kvPairsToSave = None))
      else {
        val failureStatus = handleExecutionFailure(status, None)
        retryElseFail(failureStatus)
      }
    }

  }

  // get the exit code of the job.
  def JobExitCode: Future[String] = {

    // read if the file exists
    def readRCFile(fileExists: Boolean): Future[String] =
      if (fileExists)
        asyncIo.contentAsStringAsync(jobPaths.returnCode, None, failOnOverflow = false)
      else {
        jobLogger.warn("RC file not found in aws version. Setting job to failed.")
        Future("1")
      }
    // finally : assign the yielded variable
    for {
      fileExists <- asyncIo.existsAsync(jobPaths.returnCode)
      jobRC <- readRCFile(fileExists)
    } yield jobRC
  }

  // new OOM detection
  def memoryRetryRC(job: StandardAsyncJob): Future[Boolean] = Future {
    // STATUS LOGIC:
    //   - success : container exit code is zero
    //   - command failure: container exit code > 0, no statusReason in container
    //   - OOM kill : container exit code > 0, statusReason contains "OutOfMemory" OR exit code == 137
    //   - spot kill : no container exit code set. statusReason of ATTEMPT (not container) says "host EC2 (...) terminated"
    Log.debug(s"Looking for memoryRetry in job '${job.jobId}'")
    val describeJobsResponse = batchClient.describeJobs(DescribeJobsRequest.builder.jobs(job.jobId).build)
    val jobDetail = describeJobsResponse.jobs.get(0)
    val nrAttempts = jobDetail.attempts.size
    // if job is terminated/cancelled before starting, there are no attempts.
    val lastAttemptOpt = Try(jobDetail.attempts.get(nrAttempts - 1)).toOption
    lastAttemptOpt
      .map { lastAttempt =>
        // if missing, set to failed.
        val containerRC = Try(lastAttempt.container.exitCode).toOption.getOrElse(1)

        // if not zero => get reason, else set retry to false.
        containerRC.toString match {
          case "0" =>
            Log.debug("container exit code was zero. job succeeded")
            false
          case "137" =>
            Log.info("Job failed with Container status reason : 'OutOfMemory' (code:137)")
            true
          case _ =>
            // failed job due to command errors (~ user errors) don't have a container exit reason.
            val containerStatusReason: String = {
              // if no attempts were made (rare) : container is null:
              val lastReason = Try(lastAttempt.container.reason).toOption
              if (lastReason.isEmpty) {
                log.debug("No exit reason found for container.")
              } else {
                Log.warn(s"Job failed with Container status reason : '${lastReason}'")
              }
              lastReason.getOrElse("")
            }
            // check the list of OOM-keys against the exit reason.
            val RetryMemoryKeys = memoryRetryErrorKeys.toList.flatten
            val retry = RetryMemoryKeys.exists(containerStatusReason.contains)
            Log.debug(s"Retry job based on provided keys : '${retry}'")
            retry
        }
      }
      .getOrElse {
        Log.info(s"No attempts were made for job '${job.jobId}'. no memory-related retry needed.")
        false
      }
  }

  // Despite being a "runtime" exception, BatchExceptions for 429 (too many requests) are *not* fatal:
  override def isFatal(throwable: Throwable): Boolean = throwable match {
    case be: BatchException => !be.getMessage.contains("Status Code: 429")
    case _ => super.isFatal(throwable)
  }

  override lazy val startMetadataKeyValues: Map[String, Any] =
    super[AwsBatchJobCachingActorHelper].startMetadataKeyValues

  // opportunity to send custom metadata when the run is in a terminal state, currently related to cloudwatch info
  def getTerminalMetadata(runStatus: RunStatus, jobHandle: StandardAsyncPendingExecutionHandle): Map[String, Any] = {
    // job details
    val jobDetail =
      batchClient.describeJobs(DescribeJobsRequest.builder.jobs(jobHandle.pendingJob.jobId).build).jobs.get(0)
    runStatus match {
      case _: TerminalRunStatus =>
        Map(
          AwsBatchMetadataKeys.LogStreamName -> Option(jobDetail.container.logStreamName).getOrElse("unknown"),
          AwsBatchMetadataKeys.LogGroupName -> Option(jobDetail.container.logConfiguration.options.get("awslogs-group"))
            .getOrElse("unknown"),
          // AwsBatchMetadataKeys.JobStatusReason -> Option(jobDetail.statusReason).getOrElse("unknown"),
          // region is taken by splitting the ARN of the job queue (logging is always in same region I think)
          AwsBatchMetadataKeys.LogStreamRegion -> Option(jobDetail.jobQueue.split(":")(3)).getOrElse("unknown")
        )
      case unknown => throw new RuntimeException(s"Attempt to get terminal metadata from non terminal status: $unknown")
    }
  }

  def hostAbsoluteFilePath(jobPaths: JobPaths, pathString: String): Path = {

    val pathBuilders: List[PathBuilder] = List(DefaultPathBuilder)
    val path = PathFactory.buildPath(pathString, pathBuilders)
    if (!path.isAbsolute)
      jobPaths.callExecutionRoot.resolve(path).toAbsolutePath
    else if (jobPaths.isInExecution(path.pathAsString))
      jobPaths.hostPathFromContainerPath(path.pathAsString)
    else
      jobPaths.hostPathFromContainerInputs(path.pathAsString)
  }

  override def mapOutputWomFile(womFile: WomFile): WomFile = {
    val wfile = configuration.fileSystem match {
      case AWSBatchStorageSystems.s3 =>
        womFile
      case _ =>
        val hostPath = hostAbsoluteFilePath(jobPaths, womFile.valueString)
        if (!hostPath.exists)
          throw new FileNotFoundException(s"Could not process output, file not found: ${hostPath.pathAsString}")
        womFile mapFile { _ => hostPath.pathAsString }
    }
    womFileToPath(generateAwsBatchOutputs(jobDescriptor))(wfile)
  }

  private[aws] def womFileToPath(outputs: Set[AwsBatchFileOutput])(womFile: WomFile): WomFile =
    womFile mapFile { path =>
      outputs collectFirst {
        case output if output.name == makeSafeAwsBatchReferenceName(path) => output.s3key
      } getOrElse path
    }

  override def getTerminalEvents(runStatus: RunStatus): Seq[ExecutionEvent] =
    runStatus match {
      case successStatus: RunStatus.Succeeded => successStatus.eventList
      case unknown =>
        throw new RuntimeException(s"handleExecutionSuccess not called with RunStatus.Success. Instead got $unknown")
    }

  override def retryEvaluateOutputs(exception: Exception): Boolean =
    exception match {
      case aggregated: CromwellAggregatedException =>
        aggregated.throwables.collectFirst { case s: SocketTimeoutException => s }.isDefined
      case _ => false
    }

  override def mapCommandLineWomFile(womFile: WomFile): WomFile =
    womFile.mapFile(value =>
      getPath(value) match {
        case Success(path: S3Path) => workingDisk.mountPoint.resolve(path.pathWithoutScheme).pathAsString
        case _ => value
      }
    )

  override def handlePollSuccess(oldHandle: StandardAsyncPendingExecutionHandle,
                                 state: StandardAsyncRunState
  ): Future[ExecutionHandle] = {
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
        val metadata = getTerminalMetadata(state, oldHandle)
        tellMetadata(metadata)
        handleExecutionResult(state, oldHandle)
      case s =>
        Future.successful(
          oldHandle.copy(previousState = Option(s))
        ) // Copy the current handle with updated previous status.
    }
  }

  override def handleExecutionSuccess(runStatus: StandardAsyncRunState,
                                      handle: StandardAsyncPendingExecutionHandle,
                                      returnCode: Int
  )(implicit ec: ExecutionContext): Future[ExecutionHandle] =
    evaluateOutputs() map {
      case ValidJobOutputs(outputs) =>
        // Need to make sure the paths are up to date before sending the detritus back in the response
        updateJobPaths()
        // If instance is terminated while copying stdout/stderr : status is failed while jobs outputs are ok
        //   => Retryable
        if (runStatus.toString().equals("Failed")) {
          jobLogger.warn("Got Failed RunStatus for success Execution")

          val exception = new MessageAggregation {
            override def exceptionContext: String = "Got Failed RunStatus for success Execution"
            override def errorMessages: Iterable[String] = Array("Got Failed RunStatus for success Execution")
          }
          FailedNonRetryableExecutionHandle(exception, kvPairsToSave = None)
        } else {
          SuccessfulExecutionHandle(outputs, returnCode, jobPaths.detritusPaths, getTerminalEvents(runStatus))
        }
      case InvalidJobOutputs(errors) =>
        val exception = new MessageAggregation {
          override def exceptionContext: String = "Failed to evaluate job outputs"
          override def errorMessages: Iterable[String] = errors.toList
        }
        FailedNonRetryableExecutionHandle(exception, kvPairsToSave = None)
      case JobOutputsEvaluationException(exception: Exception) if retryEvaluateOutputsAggregated(exception) =>
        // Return the execution handle in this case to retry the operation
        handle
      case JobOutputsEvaluationException(ex) => FailedNonRetryableExecutionHandle(ex, kvPairsToSave = None)
    }

  // overrides for globbing
  /**
    * Returns the shell scripting for linking a glob results file.
    *
    * @param globFile The glob.
    * @return The shell scripting.
    */
  override def globScript(globFile: WomGlobFile): String = {
    val (globName, globbedDir, _, _, _) = generateGlobPaths(globFile)
    val controlFileName = "cromwell_glob_control_file"
    val absoluteGlobValue = commandDirectory.resolve(globFile.value).pathAsString
    val globDirectory = globbedDir + "/." + globName + "/"
    val globList = globbedDir + "/." + globName + ".list"
    val globLinkCommand: String =
      (if (configuration.globLinkCommand.isDefined) {
         "( " + configuration.globLinkCommand.getOrElse("").toString + " )"

       } else {
         "( ln -L GLOB_PATTERN GLOB_DIRECTORY 2> /dev/null ) || ( ln GLOB_PATTERN GLOB_DIRECTORY )"
       }).toString
        .replaceAll("GLOB_PATTERN", absoluteGlobValue)
        .replaceAll("GLOB_DIRECTORY", globDirectory)
    // if on EFS : remove the globbing dir first, to remove leftover links from previous globs.
    val mkDirCmd: String =
      if (
        configuration.efsMntPoint.isDefined && globDirectory.startsWith(
          configuration.efsMntPoint.getOrElse("")
        )
      ) {
        jobLogger.warn("Globbing on EFS has risks.")
        jobLogger.warn(
          s"The globbing target (${globbedDir}/.${globName}/) will be overwritten when existing!"
        )
        jobLogger.warn(
          "Consider keeping globbed outputs in the cromwell-root folder"
        )
        s"rm -Rf $globDirectory $globList && mkdir"
      } else {
        "mkdir"
      }

    val controlFileContent =
      """This file is used by Cromwell to allow for globs that would not match any file.
        |By its presence it works around the limitation of some backends that do not allow empty globs.
        |Regardless of the outcome of the glob, this file will not be part of the final list of globbed files.
      """.stripMargin

    s"""|# make the directory which will keep the matching files
        |$mkDirCmd $globDirectory
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
}
