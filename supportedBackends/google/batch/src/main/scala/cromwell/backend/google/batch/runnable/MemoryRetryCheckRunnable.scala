package cromwell.backend.google.batch.runnable

import com.google.cloud.batch.v1.{Runnable, Volume}
import cromwell.backend.google.batch.api.GcpBatchRequestFactory.CreateBatchJobParameters

trait MemoryRetryCheckRunnable {

  def checkForMemoryRetryRunnables(createParameters: CreateBatchJobParameters, volumes: List[Volume]): List[Runnable] =
    createParameters.retryWithMoreMemoryKeys match {
      case Some(keys) => List(RunnableBuilder.checkForMemoryRetryRunnable(keys, volumes)).map(_.build)
      case None => List.empty[Runnable]
    }
}
