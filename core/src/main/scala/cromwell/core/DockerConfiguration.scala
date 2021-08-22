package cromwell.core

import cats.data.Validated._
import cats.syntax.apply._
import com.typesafe.config.{Config, ConfigFactory}
import common.exception.AggregatedMessageException
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import net.ceedubs.ficus.Ficus._
import org.slf4j.LoggerFactory

import scala.concurrent.duration.FiniteDuration

object DockerConfiguration {
  private val logger = LoggerFactory.getLogger("DockerConfiguration")
  lazy val dockerConfig: Config = ConfigFactory.load().getConfig("docker")
  lazy val dockerHashLookupConfig: Config = dockerConfig.getConfig("hash-lookup")

  lazy val instance: DockerConfiguration = {
    if (dockerHashLookupConfig.hasPath("gcr-api-queries-per-100-seconds")) {
      logger.warn("'docker.hash-lookup.gcr-api-queries-per-100-seconds' is no longer supported, use 'docker.hash-lookup.google.throttle' instead (see reference.conf)")
    }
    val enabled = validate { dockerHashLookupConfig.as[Boolean]("enabled") }
    val cacheEntryTtl = validate { dockerHashLookupConfig.as[FiniteDuration]("cache-entry-ttl") }
    val cacheSize = validate { dockerHashLookupConfig.as[Long]("cache-size") }
    val method: ErrorOr[DockerHashLookupMethod] = validate { dockerHashLookupConfig.as[String]("method") } map {
      case "local" => DockerLocalLookup
      case "remote" => DockerRemoteLookup
      case other => throw new IllegalArgumentException(s"Unrecognized docker hash lookup method: $other")
    }
    val sizeCompressionFactor = validate { dockerHashLookupConfig.as[Double]("size-compression-factor") }
    val maxTimeBetweenRetries = validate { dockerHashLookupConfig.as[FiniteDuration]("max-time-between-retries") }
    val maxRetries = validate { dockerHashLookupConfig.as[Int]("max-retries") }

    val dockerConfiguration = (enabled,
      cacheEntryTtl, cacheSize, method,
      sizeCompressionFactor, maxTimeBetweenRetries, maxRetries) mapN DockerConfiguration.apply

    dockerConfiguration match {
      case Valid(conf) => conf
      case Invalid(errors) => throw AggregatedMessageException("Invalid docker configuration", errors.toList)
    }
  }
}

case class DockerConfiguration(
                                enabled: Boolean,
                                cacheEntryTtl: FiniteDuration,
                                cacheSize: Long,
                                method: DockerHashLookupMethod,
                                sizeCompressionFactor: Double,
                                maxTimeBetweenRetries: FiniteDuration,
                                maxRetries: Int
                              )

sealed trait DockerHashLookupMethod

case object DockerLocalLookup extends DockerHashLookupMethod
case object DockerRemoteLookup extends DockerHashLookupMethod
