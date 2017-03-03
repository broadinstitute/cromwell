package cromwell.backend.sfs

import akka.actor.Actor
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.standard.StandardCachingActorHelper
import cromwell.core.logging.JobLogging
import cromwell.core.path.PathBuilder
import net.ceedubs.ficus.Ficus._

trait SharedFileSystemJobCachingActorHelper extends StandardCachingActorHelper {
  this: Actor with JobLogging =>

  lazy val sharedFileSystem = new SharedFileSystem {
    override lazy val pathBuilders: List[PathBuilder] = standardInitializationData.workflowPaths.pathBuilders

    override lazy val sharedFileSystemConfig: Config = {
      configurationDescriptor.backendConfig.as[Option[Config]]("filesystems.local").getOrElse(ConfigFactory.empty())
    }
  }
}
