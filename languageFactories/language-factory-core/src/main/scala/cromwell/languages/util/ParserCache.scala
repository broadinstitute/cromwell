package cromwell.languages.util

import java.util.concurrent.Callable

import cats.data.Validated.{Invalid, Valid}
import cats.instances.list._
import cats.syntax.traverse._
import com.google.common.cache.{Cache, CacheBuilder}
import com.typesafe.config.Config
import com.typesafe.scalalogging.StrictLogging
import common.validation.ErrorOr.ErrorOr
import cromwell.core.CacheConfig
import cromwell.languages.LanguageFactory
import cromwell.languages.util.ImportResolver.ImportResolver
import cromwell.languages.util.ParserCache.ParserCacheInputs
import mouse.all._
import net.ceedubs.ficus.Ficus._
import wom.core.{WorkflowSource, WorkflowUrl}
import wom.values._

import scala.concurrent.duration._

trait ParserCache[A] extends StrictLogging { this: LanguageFactory =>

  def retrieveOrCalculate(cacheInputs: ParserCacheInputs,
                          calculationCallable: Callable[ErrorOr[A]]): ErrorOr[A] = {

    (cache map { c: Cache[String, ErrorOr[A]] =>
      workflowHashKey(cacheInputs.workflowSource, cacheInputs.workflowUrl, cacheInputs.workflowRoot, cacheInputs.importResolvers) match {
        case Valid(hashKey) => c.get(hashKey, calculationCallable)
        case Invalid(errors) =>
          logger.info(s"Failed to calculate hash key for 'workflow source to WOM' cache: {}", errors.toList.mkString(", "))
          calculationCallable.call
      }
    }).getOrElse(calculationCallable.call())
  }

  private[this] def workflowHashKey(workflowSource: Option[WorkflowSource],
                                    workflowUrl: Option[WorkflowUrl],
                                    workflowRoot: Option[String],
                                    importResolvers: List[ImportResolver]): ErrorOr[String] = {
    def stringOptionToHash(opt: Option[String]): String = opt map { _.md5Sum } getOrElse ""

    val importResolversToHash: ErrorOr[String] = importResolvers.traverse(_.hashKey).map(_.mkString(":"))

    importResolversToHash map { importHash =>
      s"${stringOptionToHash(workflowSource)}:${stringOptionToHash(workflowUrl)}:${stringOptionToHash(workflowRoot)}:$importHash"
    }
  }

  private[this] lazy val cacheConfig: Option[CacheConfig] = {
    // Caching is an opt-in activity:
    for {
      _ <- enabled.option(())
      cachingConfigSection <- config.as[Option[Config]]("caching")
      cc <- CacheConfig.optionalConfig(cachingConfigSection, defaultConcurrency = 2, defaultSize = 1000L, defaultTtl = 20.minutes)
    } yield cc
  }

  private[this] lazy val cache: Option[Cache[String, ErrorOr[A]]] = cacheConfig map { c =>
    CacheBuilder.newBuilder()
      .concurrencyLevel(c.concurrency)
      .expireAfterAccess(c.ttl.length, c.ttl.unit)
      .maximumSize(c.size)
      .build[WorkflowSource, ErrorOr[A]]()
  }
}

object ParserCache {
  final case class ParserCacheInputs(workflowSource: Option[WorkflowSource],
                                     workflowUrl: Option[WorkflowUrl],
                                     workflowRoot: Option[String],
                                     importResolvers: List[ImportResolver])
}
