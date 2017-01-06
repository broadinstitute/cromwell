package cromwell.backend.sfs

import akka.actor.Actor
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.standard.{StandardCachingActorHelper, StandardInitializationData}
import cromwell.core.logging.JobLogging
import cromwell.core.path.PathBuilder
import net.ceedubs.ficus.Ficus._

trait SharedFileSystemJobCachingActorHelper extends StandardCachingActorHelper {
  this: Actor with JobLogging =>

  lazy val sharedFileSystem = new SharedFileSystem {
    override val pathBuilders: List[PathBuilder] = {
      StandardInitializationData.pathBuilders(backendInitializationDataOption)
    }
    override lazy val sharedFileSystemConfig: Config = {
      configurationDescriptor.backendConfig.as[Option[Config]]("filesystems.local").getOrElse(ConfigFactory.empty())
    }
  }
}
