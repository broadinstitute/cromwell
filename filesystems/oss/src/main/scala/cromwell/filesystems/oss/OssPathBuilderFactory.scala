package cromwell.filesystems.oss

import akka.actor.ActorSystem
import cats.syntax.apply._
import com.typesafe.config.Config
import common.validation.Validation._
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import net.ceedubs.ficus.Ficus._

import scala.concurrent.{ExecutionContext, Future}

final case class OssPathBuilderFactory(globalConfig: Config, instanceConfig: Config) extends PathBuilderFactory {
  val (endpoint, accessId, accessKey, securityToken) = (
    validate { instanceConfig.as[String]("auth.endpoint") },
    validate { instanceConfig.as[String]("auth.access-id") },
    validate { instanceConfig.as[String]("auth.access-key") },
    validate { instanceConfig.as[Option[String]]("auth.security-token") }
  ).tupled.unsafe("OSS filesystem configuration is invalid")
  
  def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext) = {
    Future.successful(OssPathBuilder.fromConfiguration(endpoint, accessId, accessKey, securityToken, options))
  }
}
