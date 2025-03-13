package cromwell.backend.impl.aws

import cromwell.backend.ReadLikeFunctions
//import cromwell.core.io.AsyncIoFunctions
//import cromwell.core.path.PathFactory
import scala.concurrent.Future
import scala.util.Try
import cromwell.core.path.{DefaultPathBuilder}
//import cromwell.backend.impl.aws.AwsBatchBackendInitializationData
//import cromwell.backend.BackendInitializationData
import cromwell.backend.standard.StandardExpressionFunctionsParams

trait AwsReadLikeFunctions extends ReadLikeFunctions {
  // standardParams for expression does not contain backend info... 
  def standardParams: StandardExpressionFunctionsParams
  

  //def backendConfig: AwsBatchConfiguration

  //val aws_config = BackendInitializationData.as[AwsBatchBackendInitializationData](standardParams.backendInitializationDataOption).configuration

  override def readFile(path: String, maxBytes: Option[Int], failOnOverflow: Boolean): Future[String] = {
    // similar to aws globbing functions, no access to the backend config is available here.... 
    ///  => using hard coded /mnt/efs path.  
    val awsPath = if (path.startsWith("/mnt/efs/")) {
        DefaultPathBuilder.get(path)
    } else {
        buildPath(path)
    }    
    Future.fromTry(Try(awsPath)) flatMap { p => asyncIo.contentAsStringAsync(p, maxBytes, failOnOverflow) }
  }
}