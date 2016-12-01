package cromwell.backend.impl.jes

import java.net.URL

import cats.data._
import cats.data.Validated._
import cats.syntax.cartesian._
import com.typesafe.config.Config
import cromwell.backend.impl.jes.authentication.JesAuths
import cromwell.core.ErrorOr._
import cromwell.filesystems.gcs.GoogleConfiguration
import lenthall.config.ValidatedConfig._
import net.ceedubs.ficus.Ficus._
import wdl4s.ExceptionWithErrors

case class JesAttributes(project: String,
                         computeServiceAccount: String,
                         auths: JesAuths,
                         executionBucket: String,
                         endpointUrl: URL,
                         maxPollingInterval: Int,
                         qps: Int)

object JesAttributes {

  private val jesKeys = Set(
    "project",
    "root",
    "maximum-polling-interval",
    "compute-service-account",
    "dockerhub",
    "genomics",
    "filesystems",
    "genomics.auth",
    "genomics.endpoint-url",
    "filesystems.gcs.auth"
  )

  private val context = "Jes"

  def apply(googleConfig: GoogleConfiguration, backendConfig: Config): JesAttributes = {
    backendConfig.warnNotRecognized(jesKeys, context)

    val project: ValidatedNel[String, String] = backendConfig.validateString("project")
    val executionBucket: ValidatedNel[String, String] = backendConfig.validateString("root")
    val endpointUrl: ErrorOr[URL] = backendConfig.validateURL("genomics.endpoint-url")
    val maxPollingInterval: Int = backendConfig.as[Option[Int]]("maximum-polling-interval").getOrElse(600)
    val computeServiceAccount: String = backendConfig.as[Option[String]]("genomics.compute-service-account").getOrElse("default")
    val genomicsAuthName: ErrorOr[String] = backendConfig.validateString("genomics.auth")
    val gcsFilesystemAuthName: ErrorOr[String] = backendConfig.validateString("filesystems.gcs.auth")

    // 1000 per 100s is the default API limit
    val qps = backendConfig.as[Option[Int]]("genomics-api-queries-per-100-seconds").getOrElse(1000) / 100

    (project |@| executionBucket |@| endpointUrl |@| genomicsAuthName |@| gcsFilesystemAuthName) map {
      (_, _, _, _, _)
    } flatMap { case (p, b, u, genomicsName, gcsName) =>
      (googleConfig.auth(genomicsName) |@| googleConfig.auth(gcsName)) map { case (genomicsAuth, gcsAuth) =>
        JesAttributes(p, computeServiceAccount, JesAuths(genomicsAuth, gcsAuth), b, u, maxPollingInterval, qps)
      }
    } match {
      case Valid(r) => r
      case Invalid(f) =>
        throw new IllegalArgumentException with ExceptionWithErrors {
          override val message = "Jes Configuration is not valid: Errors"
          override val errors = f
        }
    }
  }
}
