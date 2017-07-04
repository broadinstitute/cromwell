package cromwell.backend.impl.jes

import java.net.{URI, URL}

import cats.data.Validated._
import cats.syntax.cartesian._
import cats.syntax.validated._
import com.typesafe.config.{Config, ConfigValue}
import cromwell.backend.impl.jes.authentication.JesAuths
import cromwell.backend.impl.jes.callcaching.{CopyCachedOutputs, JesCacheHitDuplicationStrategy, UseOriginalCachedOutputs}
import cromwell.filesystems.gcs.GoogleConfiguration
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import lenthall.exception.MessageAggregation
import lenthall.validation.ErrorOr._
import lenthall.validation.Validation._
import net.ceedubs.ficus.Ficus._
import net.ceedubs.ficus.readers.{StringReader, ValueReader}
import org.slf4j.{Logger, LoggerFactory}

import scala.collection.JavaConverters._

case class JesAttributes(project: String,
                         computeServiceAccount: String,
                         auths: JesAuths,
                         executionBucket: String,
                         endpointUrl: URL,
                         maxPollingInterval: Int,
                         qps: Int Refined Positive,
                         duplicationStrategy: JesCacheHitDuplicationStrategy)

object JesAttributes {
  lazy val Logger = LoggerFactory.getLogger("JesAttributes") 
  
  val GenomicsApiDefaultQps = 1000

  private val jesKeys = Set(
    "project",
    "root",
    "maximum-polling-interval",
    "genomics.compute-service-account",
    "dockerhub",
    "dockerhub.account",
    "dockerhub.token",
    "genomics",
    "filesystems",
    "genomics.auth",
    "genomics.endpoint-url",
    "filesystems.gcs.auth",
    "filesystems.gcs.caching.duplication-strategy",
    "genomics-api-queries-per-100-seconds"
  )

  private val deprecatedJesKeys: Map[String, String] = Map(
    "genomics.default-zones" -> "default-runtime-attributes.zones"
  )

  private val context = "Jes"

  implicit val urlReader: ValueReader[URL] = StringReader.stringValueReader.map { URI.create(_).toURL }
  
  def apply(googleConfig: GoogleConfiguration, backendConfig: Config): JesAttributes = {
    val configKeys = backendConfig.entrySet().asScala.toSet map { entry: java.util.Map.Entry[String, ConfigValue] => entry.getKey }
    warnNotRecognized(configKeys, jesKeys, context, Logger)

    def warnDeprecated(keys: Set[String], deprecated: Map[String, String], context: String, logger: Logger) = {
      val deprecatedKeys = keys.intersect(deprecated.keySet)
      deprecatedKeys foreach { key => logger.warn(s"Found deprecated configuration key $key, replaced with ${deprecated.get(key)}") }
    }

    warnDeprecated(configKeys, deprecatedJesKeys, context, Logger)

    val project: ErrorOr[String] = validate { backendConfig.as[String]("project") }
    val executionBucket: ErrorOr[String] = validate { backendConfig.as[String]("root") }
    val endpointUrl: ErrorOr[URL] = validate { backendConfig.as[URL]("genomics.endpoint-url") }
    val maxPollingInterval: Int = backendConfig.as[Option[Int]]("maximum-polling-interval").getOrElse(600)
    val computeServiceAccount: String = backendConfig.as[Option[String]]("genomics.compute-service-account").getOrElse("default")
    val genomicsAuthName: ErrorOr[String] = validate { backendConfig.as[String]("genomics.auth") }
    val gcsFilesystemAuthName: ErrorOr[String] = validate { backendConfig.as[String]("filesystems.gcs.auth") }
    val qpsValidation = validateQps(backendConfig)
    val duplicationStrategy = validate { backendConfig.as[Option[String]]("filesystems.gcs.caching.duplication-strategy").getOrElse("copy") match {
      case "copy" => CopyCachedOutputs
      case "reference" => UseOriginalCachedOutputs
      case other => throw new IllegalArgumentException(s"Unrecognized caching duplication strategy: $other. Supported strategies are copy and reference. See reference.conf for more details.")
    } }


    (project |@| executionBucket |@| endpointUrl |@| genomicsAuthName |@| gcsFilesystemAuthName |@| qpsValidation |@| duplicationStrategy).tupled flatMap { 
      case (p, b, u, genomicsName, gcsName, qps, cachingStrategy) =>
      (googleConfig.auth(genomicsName) |@| googleConfig.auth(gcsName)) map { case (genomicsAuth, gcsAuth) =>
        JesAttributes(p, computeServiceAccount, JesAuths(genomicsAuth, gcsAuth), b, u, maxPollingInterval, qps, cachingStrategy)
      }
    } match {
      case Valid(r) => r
      case Invalid(f) =>
        throw new IllegalArgumentException with MessageAggregation {
          override val exceptionContext = "Jes Configuration is not valid: Errors"
          override val errorMessages = f.toList
        }
    }
  }

  def validateQps(config: Config): ErrorOr[Int Refined Positive] = {
    import eu.timepit.refined._

    val qp100s = config.as[Option[Int]]("genomics-api-queries-per-100-seconds").getOrElse(GenomicsApiDefaultQps)
    val qpsCandidate = qp100s / 100

    refineV[Positive](qpsCandidate) match {
      case Left(_) => s"Calculated QPS for Google Genomics API ($qpsCandidate/s) was not a positive integer (supplied value was $qp100s per 100s)".invalidNel
      case Right(refined) => refined.validNel
    }
  }
}
