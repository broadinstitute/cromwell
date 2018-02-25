package cromwell.core

import com.typesafe.config.ConfigFactory

import scala.concurrent.duration.FiniteDuration
import cats.data.Validated._
import cats.syntax.apply._
import common.exception.AggregatedMessageException
import common.validation.ErrorOr.ErrorOr
import net.ceedubs.ficus.Ficus._
import common.validation.Validation._

object DockerConfiguration {
  private lazy val dockerConfig = ConfigFactory.load().getConfig("docker")
  private lazy val dockerHashLookupConfig = dockerConfig.getConfig("hash-lookup")
  
  lazy val instance: DockerConfiguration = {
    val enabled = validate { dockerHashLookupConfig.as[Boolean]("enabled") }
    val gcrApiQueriesPer100Seconds = validate { dockerHashLookupConfig.as[Int]("gcr-api-queries-per-100-seconds") }
    val cacheEntryTtl = validate { dockerHashLookupConfig.as[FiniteDuration]("cache-entry-ttl") }
    val cacheSize = validate { dockerHashLookupConfig.as[Long]("cache-size") }
    val method: ErrorOr[DockerHashLookupMethod] = validate { dockerHashLookupConfig.as[String]("method") } map {
      case "local" => DockerLocalLookup
      case "remote" => DockerRemoteLookup
      case other => throw new IllegalArgumentException(s"Unrecognized docker hash lookup method: $other")
    }

    val dockerConfiguration = (enabled, gcrApiQueriesPer100Seconds, cacheEntryTtl, cacheSize, method) mapN  DockerConfiguration.apply
    
    dockerConfiguration match {
      case Valid(conf) => conf
      case Invalid(errors) => throw AggregatedMessageException("Invalid docker configuration", errors.toList)
    }
  }
}

case class DockerConfiguration(
                                enabled: Boolean,
                                gcrApiQueriesPer100Seconds: Int,
                                cacheEntryTtl: FiniteDuration,
                                cacheSize: Long,
                                method: DockerHashLookupMethod
                              )

sealed trait DockerHashLookupMethod

case object DockerLocalLookup extends DockerHashLookupMethod
case object DockerRemoteLookup extends DockerHashLookupMethod
