package cromwell.backend.impl.jes.callcaching

import cromwell.backend.BackendCacheHitCopyingActor.CopyOutputsCommand
import cromwell.backend.BackendJobExecutionActor.{JobFailedNonRetryableResponse, JobSucceededResponse}
import cromwell.backend.io.JobPaths
import cromwell.backend.standard.callcaching.{StandardCacheHitCopyingActor, StandardCacheHitCopyingActorParams}
import cromwell.core.simpleton.WdlValueBuilder
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import lenthall.util.TryUtil

import scala.util.{Failure, Success}

class JesBackendCacheHitNoCopyActor(standardParams: StandardCacheHitCopyingActorParams) extends StandardCacheHitCopyingActor(standardParams) with GcsBatchCommandBuilder {
  override def receive = {
    case CopyOutputsCommand(simpletons, jobDetritus, returnCode) => 
      val jobOutputs = WdlValueBuilder.toJobOutputs(standardParams.jobDescriptor.call.task.outputs, simpletons)
      val detritusPaths = jobDetritus.filterNot(_._1 == JobPaths.CallRootPathKey).mapValues(getPath) + (JobPaths.CallRootPathKey -> Success(jobPaths.callRoot))
      
      val response = TryUtil.sequenceMap(detritusPaths) match {
        case Success(detritus) => JobSucceededResponse(standardParams.jobDescriptor.key, returnCode, jobOutputs, Option(detritus), Seq.empty, None)
        case Failure(failure) => JobFailedNonRetryableResponse(jobDescriptor.key, failure, None)
      }
      
      sender() ! response
      context stop self
  }
}
