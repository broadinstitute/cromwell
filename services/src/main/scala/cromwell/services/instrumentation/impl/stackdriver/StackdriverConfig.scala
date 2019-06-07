package cromwell.services.instrumentation.impl.stackdriver

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import com.typesafe.config.Config
import common.exception.MessageAggregation
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation.validate
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.{GoogleAuthMode, UserServiceAccountMode}
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._

case class StackdriverConfig(googleProject: String,
                             auth: GoogleAuthMode,
                             flushRate: FiniteDuration,
                             cromwellInstanceIdentifier: Option[String],
                             cromwellInstanceRole: Option[String],
                             cromwellPerfTestCase: Option[String])

object StackdriverConfig {
  val CromwellInstanceIdentifier = "cromwell-instance-identifier"
  val CromwellInstanceRole = "cromwell-instance-role"
  val CromwellPerfTest = "cromwell-perf-test-case"

  private def validateFlushRate(rate: FiniteDuration): ErrorOr[FiniteDuration] = {
    validate[FiniteDuration](rate) match {
      case Valid(duration) => duration match {
        case d: FiniteDuration if d.<(1.minute) => (s"`flush-rate` can't be less than 1 minute. Current rate is `${d.toString}`. " +
          s"Google's Stackdriver API needs each metric to be sent not more than once per minute.").invalidNel
        case d => d.validNel
      }
      case Invalid(e) => e.invalid
    }
  }


  private def validateAuth(authSchemeString: String, googleConfiguration: GoogleConfiguration): ErrorOr[GoogleAuthMode] = {
    validate[String](authSchemeString) match {
      case Valid(schemeString) => googleConfiguration.auth(schemeString) match {
        case Valid(_: UserServiceAccountMode) => s"`auth` scheme: $schemeString is of type `user_service_account` which is not allowed for Stackdriver instrumentation.".invalidNel
        case Valid(auth) => auth.valid
        case Invalid(error) => s"`auth` scheme is invalid. Errors: ${error.toString}".invalidNel
      }
      case Invalid(e) => e.invalid
    }
  }


  def apply(serviceConfig: Config, globalConfig: Config): StackdriverConfig = {
    val googleConfiguration: GoogleConfiguration = GoogleConfiguration(globalConfig)

    val googleProject: ErrorOr[String] = validate[String] { serviceConfig.as[String]("google-project") }
    val authScheme: ErrorOr[GoogleAuthMode] = validateAuth(serviceConfig.as[String]("auth"), googleConfiguration)
    val flushRate: ErrorOr[FiniteDuration] = validateFlushRate(serviceConfig.as[FiniteDuration]("flush-rate"))
    val cromwellInstanceIdentifier: ErrorOr[Option[String]] = serviceConfig.as[Option[String]](CromwellInstanceIdentifier).validNel
    val cromwellInstanceRole: ErrorOr[Option[String]] = serviceConfig.as[Option[String]](CromwellInstanceRole).validNel
    val cromwellPerfTestCase: ErrorOr[Option[String]] = serviceConfig.as[Option[String]](CromwellPerfTest).validNel

    (googleProject, authScheme, flushRate, cromwellInstanceIdentifier, cromwellInstanceRole, cromwellPerfTestCase).mapN({ (p, a, f, i, r, t) =>
      new StackdriverConfig(p, a, f, i, r, t)
    }).valueOr(errors => throw new IllegalArgumentException with MessageAggregation {
      override val exceptionContext = "Stackdriver instrumentation config is invalid."
      override val errorMessages = errors.toList
    })
  }
}
