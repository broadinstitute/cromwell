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

import java.net.SocketTimeoutException

import cats.data.Validated.Valid
import common.util.StringUtil._
import common.validation.Validation._
import cromwell.backend._
import cromwell.backend.async.{ExecutionHandle, PendingExecutionHandle}
import cromwell.backend.impl.aws.RunStatus.TerminalRunStatus
import cromwell.backend.impl.aws.io._
import cromwell.backend.io.DirectoryFunctions
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core._
import cromwell.core.logging.JobLogger
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.filesystems.s3.S3Path
import cromwell.filesystems.s3.batch.S3BatchCommandBuilder
import cromwell.services.keyvalue.KeyValueServiceActor._
import cromwell.services.keyvalue.KvClient
import org.slf4j.LoggerFactory
import wom.CommandSetupSideEffectFile
import wom.callable.Callable.OutputDefinition
import wom.core.FullyQualifiedName
import wom.expression.NoIoFunctionSet
import wom.types.{WomArrayType, WomSingleFileType}
import wom.values._

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Success, Try}

object AwsBatchAsyncBackendJobExecutionActor {
  val AwsBatchOperationIdKey = "__aws_batch_operation_id"

  type AwsBatchPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, AwsBatchJob, RunStatus]
}

class AwsBatchAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
  extends BackendJobLifecycleActor with StandardAsyncExecutionActor with AwsBatchJobCachingActorHelper
    with KvClient {

  override lazy val ioCommandBuilder = S3BatchCommandBuilder

  import AwsBatchAsyncBackendJobExecutionActor._

  val Log = LoggerFactory.getLogger(AwsBatchAsyncBackendJobExecutionActor.getClass)
  val logger = new JobLogger("AwsBatchRun", jobDescriptor.workflowDescriptor.id, jobDescriptor.key.tag, None, Set(Log))

  override type StandardAsyncRunInfo = AwsBatchJob

  override type StandardAsyncRunStatus = RunStatus

  override lazy val pollBackOff = SimpleExponentialBackoff(1.second, 5.minutes, 1.1)

  override lazy val executeOrRecoverBackOff = SimpleExponentialBackoff(
    initialInterval = 3 seconds, maxInterval = 20 seconds, multiplier = 1.1)

  private lazy val jobDockerImage = jobDescriptor.maybeCallCachingEligible.dockerHash.getOrElse(runtimeAttributes.dockerImage)

  override lazy val dockerImageUsed: Option[String] = Option(jobDockerImage)

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
  lazy val batchJob = {
    val script = commandScriptContents match {
      case Valid(s) => s
      case errors => new RuntimeException(errors.toList.mkString(", "))
    }
    AwsBatchJob(jobDescriptor, runtimeAttributes,
                                  instantiatedCommand.commandString,
                                  script.toString,
                                  rcPath.toString, executionStdout, executionStderr,
                                  Seq.empty[AwsBatchParameter])
  }
  /* Tries to abort the job in flight
   *
   * @param job A StandardAsyncJob object (has jobId value) to cancel
   * @return Nothing
   *
   */
  override def tryAbort(job: StandardAsyncJob): Unit = {
    batchJob.abort(job.jobId) // job.JobId should be the AWS Batch Job Id based on analysis of other backends
    Log.info(s"Attempted CancelJob operation in AWS Batch for Job ID ${job.jobId}. There were no errors during the operation")
    Log.info(s"We have normality. Anything you still can't cope with is therefore your own problem")
    Log.info(s"https://www.youtube.com/watch?v=YCRxnjE7JVs")
    ()
  }

  override def requestsAbortAndDiesImmediately: Boolean = false

  /**
    * Takes two arrays of remote and local WOM File paths and generates the necessary AwsBatchInputs.
    */
  private def inputsFromWomFiles(namePrefix: String,
                                    remotePathArray: Seq[WomFile],
                                    localPathArray: Seq[WomFile],
                                    jobDescriptor: BackendJobDescriptor): Iterable[AwsBatchInput] = {
    (remotePathArray zip localPathArray zipWithIndex) flatMap {
      case ((remotePath, localPath), index) =>
        Seq(AwsBatchFileInput(s"$namePrefix-$index", remotePath.valueString, DefaultPathBuilder.get(localPath.valueString), workingDisk))
    }
  }

  /**
    * Turns WomFiles into relative paths.  These paths are relative to the working disk.
    *
    * relativeLocalizationPath("foo/bar.txt") -> "foo/bar.txt"
    * relativeLocalizationPath("s3://some/bucket/foo.txt") -> "some/bucket/foo.txt"
    */
  private def relativeLocalizationPath(file: WomFile): WomFile = {
    file.mapFile(value =>
      getPath(value) match {
        case Success(path) => path.pathWithoutScheme
        case _ => value
      }
    )
  }

  private[aws] def generateAwsBatchInputs(jobDescriptor: BackendJobDescriptor): Set[AwsBatchInput] = {
    // We need to tell PAPI about files that were created as part of command instantiation (these need to be defined
    // as inputs that will be localized down to the VM). Make up 'names' for these files that are just the short
    // md5's of their paths.
    val writeFunctionFiles = instantiatedCommand.createdFiles map { f => f.file.value.md5SumShort -> List(f) } toMap

    def localizationPath(f: CommandSetupSideEffectFile) =
      f.relativeLocalPath.fold(ifEmpty = relativeLocalizationPath(f.file))(WomFile(f.file.womFileType, _))
    val writeFunctionInputs = writeFunctionFiles flatMap {
      case (name, files) => inputsFromWomFiles(name, files.map(_.file), files.map(localizationPath), jobDescriptor)
    }

    // Collect all WomFiles from inputs to the call.
    val callInputFiles: Map[FullyQualifiedName, Seq[WomFile]] = jobDescriptor.fullyQualifiedInputs mapValues {
      womFile =>
        val arrays: Seq[WomArray] = womFile collectAsSeq {
          case womFile: WomFile =>
            val files: List[WomSingleFile] = DirectoryFunctions
              .listWomSingleFiles(womFile, callPaths.workflowPaths)
              .toTry(s"Error getting single files for $womFile").get
            WomArray(WomArrayType(WomSingleFileType), files)
        }

        arrays.flatMap(_.value).collect {
          case womFile: WomFile => womFile
        }
    } map identity // <-- unlazy the mapValues

    val callInputInputs = callInputFiles flatMap {
      case (name, files) => inputsFromWomFiles(name, files, files.map(relativeLocalizationPath), jobDescriptor)
    }

    (writeFunctionInputs ++ callInputInputs).toSet

  }

  /**
    * Given a path (relative or absolute), returns a (Path, AwsBatchVolume) tuple where the Path is
    * relative to the Volume's mount point
    *
    * @throws Exception if the `path` does not live in one of the supplied `disks`
    */
  private def relativePathAndVolume(path: String, disks: Seq[AwsBatchVolume]): (Path, AwsBatchVolume) = {
    val absolutePath = DefaultPathBuilder.get(path) match {
      case p if !p.isAbsolute => AwsBatchWorkingDisk.MountPoint.resolve(p)
      case p => p
    }

    disks.find(d => absolutePath.startsWith(d.mountPoint)) match {
      case Some(disk) => (disk.mountPoint.relativize(absolutePath), disk)
      case None =>
        throw new Exception(s"Absolute path $path doesn't appear to be under any mount points: ${disks.map(_.toString).mkString(", ")}")
    }
  }

  private def makeSafeAwsBatchReferenceName(referenceName: String) = {
    if (referenceName.length <= 127) referenceName else referenceName.md5Sum
  }

  private[aws] def generateAwsBatchOutputs(jobDescriptor: BackendJobDescriptor): Set[AwsBatchFileOutput] = {
    import cats.syntax.validated._
    def evaluateFiles(output: OutputDefinition): List[WomFile] = {
      Try(
        output.expression.evaluateFiles(jobDescriptor.localInputs, NoIoFunctionSet, output.womType).map(_.toList map { _.file })
      ).getOrElse(List.empty[WomFile].validNel)
       .getOrElse(List.empty)
    }

    val womFileOutputs = jobDescriptor.taskCall.callable.outputs.flatMap(evaluateFiles) map relativeLocalizationPath

    val outputs: Seq[AwsBatchFileOutput] = womFileOutputs.distinct flatMap {
      _.flattenFiles flatMap {
        case unlistedDirectory: WomUnlistedDirectory => generateUnlistedDirectoryOutputs(unlistedDirectory)
        case singleFile: WomSingleFile => generateAwsBatchSingleFileOutputs(singleFile)
        case globFile: WomGlobFile => generateAwsBatchGlobFileOutputs(globFile)
      }
    }

    val additionalGlobOutput = jobDescriptor.taskCall.callable.additionalGlob.toList.flatMap(generateAwsBatchGlobFileOutputs).toSet

    outputs.toSet ++ additionalGlobOutput
  }

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
        directoryDisk
      ),
      // The collection list file:
      AwsBatchFileOutput(
        makeSafeAwsBatchReferenceName(directoryPath),
        dirDestinationPath,
        DefaultPathBuilder.get(directoryPath + "*"),
        directoryDisk
      )
    )
  }

  private def generateAwsBatchSingleFileOutputs(womFile: WomSingleFile): List[AwsBatchFileOutput] = {
    val destination = callRootPath.resolve(womFile.value.stripPrefix("/")).pathAsString
    val (relpath, disk) = relativePathAndVolume(womFile.value, runtimeAttributes.disks)
    val output = AwsBatchFileOutput(makeSafeAwsBatchReferenceName(womFile.value), destination, relpath, disk)
    List(output)
  }

  private def generateAwsBatchGlobFileOutputs(womFile: WomGlobFile): List[AwsBatchFileOutput] = {
    val globName = GlobFunctions.globName(womFile.value)
    val globDirectory = globName + "/"
    val globListFile = globName + ".list"
    val globDirectoryDestinationPath = callRootPath.resolve(globDirectory).pathAsString
    val globListFileDestinationPath = callRootPath.resolve(globListFile).pathAsString

    val (_, globDirectoryDisk) = relativePathAndVolume(womFile.value, runtimeAttributes.disks)

    // We need both the glob directory and the glob list:
    List(
      // The glob directory:
      AwsBatchFileOutput(makeSafeAwsBatchReferenceName(globDirectory), globDirectoryDestinationPath, DefaultPathBuilder.get(globDirectory + "*"), globDirectoryDisk),
      // The glob list file:
      AwsBatchFileOutput(makeSafeAwsBatchReferenceName(globListFile), globListFileDestinationPath, DefaultPathBuilder.get(globListFile), globDirectoryDisk)
    )
  }

  override lazy val commandDirectory: Path = AwsBatchWorkingDisk.MountPoint

  override def globParentDirectory(womGlobFile: WomGlobFile): Path = {
    val (_, disk) = relativePathAndVolume(womGlobFile.value, runtimeAttributes.disks)
    disk.mountPoint
  }

  override def isTerminal(runStatus: RunStatus): Boolean = {
    runStatus match {
      case _: TerminalRunStatus => true
      case _ => false
    }
  }

  // Primary entry point for cromwell to actually run something
  override def executeAsync(): Future[ExecutionHandle] = {
    for {
      submitJobResponse <- Future.fromTry(batchJob.submitJob())
    } yield PendingExecutionHandle(jobDescriptor, StandardAsyncJob(submitJobResponse.jobId), Option(batchJob), previousStatus = None)
  }

  val futureKvJobKey = KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt + 1)

  // This is called by Cromwell after initial execution (see executeAsync above)
  // It expects a Future[RunStatus]. In this case we'll simply call the
  // AWS Batch API to do this. The AwsBatchJob object in the PendingExecutionHandle
  // will have the actual status call, so our job here is simply to pull the
  // object out of the handle and execute the underlying status method
  override def pollStatusAsync(handle: AwsBatchPendingExecutionHandle): Future[RunStatus] = {
    val jobId = handle.pendingJob.jobId
    val job = handle.runInfo match {
      case Some(job) => job
      case None => throw new RuntimeException(s"pollStatusAsync called but job not available. This should not happen. Jobid ${jobId}")
    }
    val status = job.status(jobId)
    status collect { case _: TerminalRunStatus => writeToPaths(jobId, job) }
    Future.fromTry(status)
  }

  /* Writes job output to location for Cromwell to run.
   * This should be in S3. There is also some Cromwell lifecycle magic to handle
   * something here, so I don't think this is the right location. The
   * intent I believe is that the job writes it's own output, error, etc to
   * S3, and through the magic that is the NIO Filesystem, we just automagically
   * pick it up from there. However, AFAIK, there is no NIO Filesystem for
   * AWS SDK v2 (preview 9). So, in AwsWorkflowPaths I've disconnected the
   * S3 PathBuilder stuff in favor of just doing it like TES does, which is a
   * local filesystem path. We'll pick up job information from CloudWatch Logs
   * and the describejobs API call (via our AwsBatchJob object). This will
   * need to be undone at some point and wired in more like the other backends.
   *
   *
   */
  private def writeToPaths(jobId: String, job: AwsBatchJob) = {
    // import cromwell.filesystems.s3.S3PathBuilder
    val detail = job.detail(jobId)
    // TODO: Apparently in a bad config job.rc(detail) can be null, which
    //       ends up not writing the returnCode file. This fails the job, for
    //       a failed job, so not sure it's a real problem at this time
    Log.info("Job Complete. Exit code: " + job.rc(detail))
    Log.info("Output path: " + jobPaths.standardPaths.output.toString)
    val (rc, stdout, stderr) = AwsBatchJob.parseOutput(job.output(detail))
    // TODO: Is this the right place to be writing out to S3?
    jobPaths.standardPaths.output.write(stdout)
    jobPaths.standardPaths.error.write(stderr)
    jobPaths.returnCode.write(rc.toString)
  }


  override lazy val startMetadataKeyValues: Map[String, Any] = super[AwsBatchJobCachingActorHelper].startMetadataKeyValues

  override def getTerminalMetadata(runStatus: RunStatus): Map[String, Any] = {
    runStatus match {
      case _: TerminalRunStatus => Map()
      case unknown => throw new RuntimeException(s"Attempt to get terminal metadata from non terminal status: $unknown")
    }
  }

  override def mapOutputWomFile(womFile: WomFile): WomFile = {
    womFileToPath(generateAwsBatchOutputs(jobDescriptor))(womFile)
  }

  private[aws] def womFileToPath(outputs: Set[AwsBatchFileOutput])(womFile: WomFile): WomFile = {
    womFile mapFile { path =>
      outputs collectFirst {
        case output if output.name == makeSafeAwsBatchReferenceName(path) => output.s3key
      } getOrElse path
    }
  }

  override def getTerminalEvents(runStatus: RunStatus): Seq[ExecutionEvent] = {
    runStatus match {
      case successStatus: RunStatus.Succeeded => successStatus.eventList
      case unknown =>
        throw new RuntimeException(s"handleExecutionSuccess not called with RunStatus.Success. Instead got $unknown")
    }
  }

  override def retryEvaluateOutputs(exception: Exception): Boolean = {
    exception match {
      case aggregated: CromwellAggregatedException =>
        aggregated.throwables.collectFirst { case s: SocketTimeoutException => s }.isDefined
      case _ => false
    }
  }

  override def mapCommandLineWomFile(womFile: WomFile): WomFile = {
    womFile.mapFile(value =>
      getPath(value) match {
        case Success(path: S3Path) => workingDisk.mountPoint.resolve(path.pathWithoutScheme).pathAsString
        case _ => value
      }
    )
  }
}
