package cromwell.backend.impl.jes

import java.net.URL

import com.typesafe.config.Config
import cromwell.backend.impl.jes.JesImplicits.GoogleAuthWorkflowOptions
import cromwell.core.{ErrorOr, WorkflowOptions}
import cromwell.filesystems.gcs.{GoogleAuthMode, GoogleConfiguration}
import lenthall.config.ScalaConfig._
import lenthall.config.ValidatedConfig._
import wdl4s.ExceptionWithErrors

import scala.language.postfixOps
import scalaz.Scalaz._
import scalaz.Validation.FlatMap._
import scalaz._

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

    val project: ErrorOr[String] = backendConfig.validateString("project")
    val executionBucket: ErrorOr[String] = backendConfig.validateString("root")
    val endpointUrl: ErrorOr[URL] = backendConfig.validateURL("genomics.endpoint-url")
    val maxPollingInterval: Int = backendConfig.getIntOption("maximum-polling-interval").getOrElse(600)
    val genomicsAuthName: ErrorOr[String] = backendConfig.validateString("genomics.auth")
    val gcsFilesystemAuthName: ErrorOr[String] = backendConfig.validateString("filesystems.gcs.auth")

    (project |@| executionBucket |@| endpointUrl |@| genomicsAuthName |@| gcsFilesystemAuthName) {
      (_, _, _, _, _)
    } flatMap { case (p, b, u, genomicsName, gcsName) =>
      (googleConfig.auth(genomicsName) |@| googleConfig.auth(gcsName)) { case (genomicsAuth, gcsAuth) =>
        JesAttributes(p, genomicsAuth, gcsAuth, b, u, maxPollingInterval)
      }
    } match {
      case Success(r) => r
      case Failure(f) =>
        throw new IllegalArgumentException with ExceptionWithErrors {
          override val message = "Jes Configuration is not valid: Errors"
          override val errors = f
        }
    }
  }
}
