package cromwell.backend.impl.sfs.config

import akka.actor.Props
import com.typesafe.config.Config
import cromwell.backend.standard.callcaching.StandardFileHashingActor.SingleFileHashRequest
import cromwell.backend.standard.callcaching.{StandardFileHashingActor, StandardFileHashingActorParams}
import cromwell.filesystems.gcs.GcsBatchCommandBuilder
import net.ceedubs.ficus.Ficus._

import scala.util.Try

object ConfigBackendFileHashingActor {
  def props(standardParams: StandardFileHashingActorParams) = Props(new ConfigBackendFileHashingActor(standardParams))
}

class ConfigBackendFileHashingActor(standardParams: StandardFileHashingActorParams) extends StandardFileHashingActor(standardParams) with GcsBatchCommandBuilder {

  lazy val hashingStrategy: ConfigHashingStrategy = {
    configurationDescriptor.backendConfig.as[Option[Config]]("filesystems.local.caching") map ConfigHashingStrategy.apply getOrElse ConfigHashingStrategy.defaultStrategy
  }
  
  override def customHashStrategy(fileRequest: SingleFileHashRequest): Option[Try[String]] = {
    log.debug(hashingStrategy.toString)
    Option(hashingStrategy.getHash(fileRequest, log))
  }
}