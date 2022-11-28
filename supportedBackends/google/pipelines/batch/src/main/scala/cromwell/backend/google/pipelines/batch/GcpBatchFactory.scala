package cromwell.backend.google.pipelines.batch

import com.google.cloud.batch.v1.BatchServiceClient
//import cromwell.backend.standard.StandardAsyncJob

case class GcpBatchFactory(applicationName: String) {

  val batchServiceClient = BatchServiceClient.create

  //def getRequest(job: StandardAsyncJob)


}
