package cromwell.backend.impl.jes

import java.net.URL

import cats.data._
import cats.data.Validated._
import cats.syntax.cartesian._
import cats.syntax.validated._
import com.typesafe.config.Config
import cromwell.backend.impl.jes.JesImplicits.GoogleAuthWorkflowOptions
import cromwell.core.WorkflowOptions
import cromwell.filesystems.gcs.{GoogleAuthMode, GoogleConfiguration}
import lenthall.config.ScalaConfig._
import lenthall.config.ValidatedConfig._
import cromwell.core.ErrorOr._
import wdl4s.ExceptionWithErrors

import scala.language.postfixOps

case class JesAttributes(project: String,
                         genomicsAuth: GoogleAuthMode,
                         gcsFilesystemAuth: GoogleAuthMode,
                         executionBucket: String,
                         endpointUrl: URL,
                         maxPollingInterval: Int) {
  def genomicsCredential(options: WorkflowOptions) = genomicsAuth.credential(options.toGoogleAuthOptions)
  def gcsCredential(options: WorkflowOptions) = gcsFilesystemAuth.credential(options.toGoogleAuthOptions)
}

object JesAttributes {

  private val jesKeys = Set(
    "project",
    "root",
    "maximum-polling-interval",
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
    val maxPollingInterval: Int = backendConfig.getIntOption("maximum-polling-interval").getOrElse(600)
    val genomicsAuthName: ErrorOr[String] = backendConfig.validateString("genomics.auth")
    val gcsFilesystemAuthName: ErrorOr[String] = backendConfig.validateString("filesystems.gcs.auth")

    (project |@| executionBucket |@| endpointUrl |@| genomicsAuthName |@| gcsFilesystemAuthName) map {
      (_, _, _, _, _)
    } flatMap { case (p, b, u, genomicsName, gcsName) =>
      (googleConfig.auth(genomicsName) |@| googleConfig.auth(gcsName)) map { case (genomicsAuth, gcsAuth) =>
        JesAttributes(p, genomicsAuth, gcsAuth, b, u, maxPollingInterval)
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
