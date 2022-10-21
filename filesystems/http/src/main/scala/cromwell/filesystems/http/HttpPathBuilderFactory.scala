package cromwell.filesystems.http

import akka.actor.{ActorRef, ActorSystem}
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.{PathBuilder, PathBuilderFactory}

import scala.concurrent.{ExecutionContext, Future}

class HttpPathBuilderFactory(globalConfig: Config, instanceConfig: Config) extends PathBuilderFactory {
  override def withOptions(options: WorkflowOptions, serviceRegistryActor: ActorRef)
                          (implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder] = {
    Future {
      new HttpPathBuilder
    }
  }
}
