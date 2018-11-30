//package cromwell.filesystems.demo.dos
//
//import akka.actor.ActorSystem
//import com.typesafe.config.Config
//import cromwell.core.WorkflowOptions
//import cromwell.core.path.{PathBuilder, PathBuilderFactory}
//
//import scala.concurrent.{ExecutionContext, Future}
//
//class DemoDosPathBuilderFactory(globalConfig: Config, instanceConfig: Config) extends PathBuilderFactory {
//  override def withOptions(options: WorkflowOptions)
//                          (implicit as: ActorSystem, ec: ExecutionContext): Future[PathBuilder] = {
//    Future {
//      new DemoDosPathBuilder
//    }
//  }
//}
