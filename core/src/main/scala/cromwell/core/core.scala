package cromwell.core

import com.typesafe.config.Config
import common.exception.ThrowableAggregation
import cromwell.core.path.Path
import mouse.boolean._

import scala.concurrent.duration.FiniteDuration


case class StandardPaths(output: Path, error: Path)

case class CallContext(root: Path, standardPaths: StandardPaths, isDocker: Boolean)

/**  Marker trait for Cromwell exceptions that are to be treated as fatal (non-retryable) */
trait CromwellFatalExceptionMarker { this: Throwable => }

object CromwellFatalException {
  // Don't wrap if it's already a fatal exception
  def apply(throwable: Throwable) = throwable match {
    case e: CromwellFatalExceptionMarker => e
    case e => new CromwellFatalException(e)
  }
  def unapply(e: CromwellFatalException): Option[Throwable] = Option(e.exception)
}

class CromwellFatalException(val exception: Throwable) extends Exception(exception) with CromwellFatalExceptionMarker

case class CromwellAggregatedException(throwables: Seq[Throwable], exceptionContext: String = "")
  extends Exception with ThrowableAggregation

case class CacheConfig(concurrency: Int, size: Long, ttl: FiniteDuration)

import net.ceedubs.ficus.Ficus._
object CacheConfig {
  /**
    * From an optional `Config` entry and specified defaults, always return a `CacheConfig` object.
    */
  def config(caching: Option[Config], defaultConcurrency: Int, defaultSize: Long, defaultTtl: FiniteDuration): CacheConfig = {
    caching flatMap { c =>
      optionalConfig(c, defaultConcurrency = defaultConcurrency, defaultSize = defaultSize, defaultTtl = defaultTtl)
    } getOrElse CacheConfig(concurrency = defaultConcurrency, size = defaultSize, ttl = defaultTtl)
  }

  /**
    * From a non-optional `Config` and specified defaults, if caching is enabled return a `CacheConfig` object wrapped in a `Some`,
    * otherwise return `None`.
    */
  def optionalConfig(caching: Config, defaultConcurrency: Int, defaultSize: Long, defaultTtl: FiniteDuration): Option[CacheConfig] = {
    val cachingEnabled = caching.getOrElse("enabled", false)

    cachingEnabled.option(
      CacheConfig(
        concurrency = caching.getOrElse("concurrency", defaultConcurrency),
        size = caching.getOrElse("size", defaultSize),
        ttl = caching.getOrElse("ttl", defaultTtl)
      )
    )
  }
}
