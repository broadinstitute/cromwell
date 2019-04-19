package cromwell.core

import com.typesafe.config.ConfigFactory

import scala.concurrent.duration.FiniteDuration
import cats.data.Validated._
import cats.syntax.apply._
import common.exception.AggregatedMessageException
import common.validation.ErrorOr.ErrorOr
import net.ceedubs.ficus.Ficus._
import common.validation.Validation._
import cromwell.core.io.Throttle
import org.slf4j.LoggerFactory

import scala.concurrent.duration._

object DockerConfiguration {
  private val logger = LoggerFactory.getLogger("DockerConfiguration")
  lazy val dockerConfig = ConfigFactory.load().getConfig("docker")
  lazy val dockerHashLookupConfig = dockerConfig.getConfig("hash-lookup")

  lazy val instance: DockerConfiguration = {
    if (dockerHashLookupConfig.hasPath("gcr-api-queries-per-100-seconds")) {
      logger.warn("'docker.hash-lookup.gcr-api-queries-per-100-seconds' is being deprecated, use 'docker.hash-lookup.gcr.throttle' instead (see reference.conf)")
    }
    val enabled = validate { dockerHashLookupConfig.as[Boolean]("enabled") }
    // Deprecate but keep for backwards compatibility
    val gcrApiQueriesPer100Seconds = validate { 
      dockerHashLookupConfig.getAs[Int]("gcr-api-queries-per-100-seconds") map { v =>
        Throttle(v, 100.seconds, v)
      }
    }
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

    val dockerConfiguration = (enabled, gcrApiQueriesPer100Seconds,
      cacheEntryTtl, cacheSize, method,
      sizeCompressionFactor, maxTimeBetweenRetries, maxRetries) mapN  DockerConfiguration.apply

    dockerConfiguration match {
      case Valid(conf) => conf
      case Invalid(errors) => throw AggregatedMessageException("Invalid docker configuration", errors.toList)
    }
  }
}

case class DockerConfiguration(
                                enabled: Boolean,
                                deprecatedGcrApiQueriesPer100Seconds: Option[Throttle],
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
