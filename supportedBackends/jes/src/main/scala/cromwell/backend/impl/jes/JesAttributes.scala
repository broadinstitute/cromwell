package cromwell.backend.impl.jes

import java.net.URL

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.core.{ErrorOr, WorkflowOptions}
import cromwell.filesystems.gcs.{GoogleAuthMode, GoogleConfiguration}
import cromwell.backend.impl.jes.JesImplicits.GoogleAuthWorkflowOptions
import lenthall.config.ScalaConfig._
import lenthall.config.ValidatedConfig._
import wdl4s.ThrowableWithErrors

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

  def assertWorkflowOptions(options: WorkflowOptions): Unit = {
    // These methods throw on bad options
    genomicsAuth.assertWorkflowOptions(options.toGoogleAuthOptions)
    gcsFilesystemAuth.assertWorkflowOptions(options.toGoogleAuthOptions)
  }
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

  def apply(googleConfig: GoogleConfiguration, configurationDescriptor: BackendConfigurationDescriptor): JesAttributes = {
    val backendConfig = configurationDescriptor.backendConfig
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
        throw new IllegalArgumentException() with ThrowableWithErrors {
          override val message = "Jes Configuration is not valid: Errors"
          override val errors = f
        }
    }
  }
}
