package cromwell.filesystems.minio

import com.typesafe.config.Config
import cromwell.core.path.PathBuilderFactory
import akka.actor.ActorSystem
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilder
import scala.concurrent.{ExecutionContext, Future}


class MinioPathBuilderFactory (globalConfig: Config, instanceConfig: Config) extends PathBuilderFactory {

  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder] = {
      println("Miniopath builder: " + options)
    Future {
      MinioPathBuilder
    }
  }
  
}
