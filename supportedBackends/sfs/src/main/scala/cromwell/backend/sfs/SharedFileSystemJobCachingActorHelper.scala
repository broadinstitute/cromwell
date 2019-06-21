package cromwell.backend.sfs

import akka.actor.{Actor, ActorContext}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.standard.StandardCachingActorHelper
import cromwell.core.logging.JobLogging
import cromwell.core.path.{Path, PathBuilder}
import net.ceedubs.ficus.Ficus._

trait SharedFileSystemJobCachingActorHelper extends StandardCachingActorHelper {
  this: Actor with JobLogging =>

  lazy val sharedFileSystem = new SharedFileSystem {

    override implicit def actorContext: ActorContext = context

    override lazy val pathBuilders: List[PathBuilder] = standardInitializationData.workflowPaths.pathBuilders

    override lazy val sharedFileSystemConfig: Config = {
      configurationDescriptor.backendConfig.as[Option[Config]]("filesystems.local").getOrElse(ConfigFactory.empty())
    }

    // cachedCopyDir should be on the same physical filesystem as the execution root.
    // WDL workflow names may not contain '-' so using 'cached-inputs' will certainly
    // not collide with any workflows in the root directory.
    override lazy val cachedCopyDir: Option[Path] = Option(
      workflowPaths.executionRoot.createChild("cached-inputs", asDirectory = true))
  }
}
