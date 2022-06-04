package cromwell.filesystems.s3

import akka.actor.ActorSystem
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import net.ceedubs.ficus.Ficus._

import scala.concurrent.{ExecutionContext, Future}

// The constructor of this class is required to be Config, Config by cromwell
// So, we need to take this config and get the AuthMode out of it
final case class GenericS3PathBuilderFactory private(globalConfig: Config, instanceConfig: Config)
  extends PathBuilderFactory {

  // Grab the authMode out of configuration
  val authModeAsString: String = instanceConfig.as[String]("auth")

  val config = globalConfig.getConfig("genericS3")

  val endpointString: String = config.getString("endpoint")
  val accessString: String = config.getString("access-key")
  val secretString: String = config.getString("secret-key")

  def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[GenericS3PathBuilder] = {
    Future( GenericS3PathBuilder.fromStaticMode(endpointString, accessString, secretString) )
  }

}

object GenericS3PathBuilderFactory {
  def apply(globalConfig: Config, instanceConfig: Config): GenericS3PathBuilderFactory = {
    new GenericS3PathBuilderFactory(globalConfig, instanceConfig)
  }
}
