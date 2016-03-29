package cromwell.engine.backend.jes

import java.net.URL

import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.ErrorOr
import lenthall.config.ScalaConfig._
import lenthall.config.ValidatedConfig._
import wdl4s.ThrowableWithErrors

import scala.language.postfixOps
import scalaz.Scalaz._
import scalaz._

case class JesAttributes(project: String, executionBucket: String, endpointUrl: URL, maxPollingInterval: Int)

object JesAttributes {

  private val jesKeys = Set(
    "project",
    "root",
    "endpointUrl",
    "maximumPollingInterval"
  )

  private val context = "Jes"

  def apply(config: Config): JesAttributes = {

    config.warnNotRecognized(jesKeys, context)

    val project: ErrorOr[String] = config.validateString("project")
    val executionBucket: ErrorOr[String] = config.validateString("root")
    val endpointUrl: ErrorOr[URL] = config.validateURL("endpointUrl")
    val maxPollingInterval: Int = config.getIntOption("maximumPollingInterval").getOrElse(600)

    (project |@| executionBucket |@| endpointUrl) {
      JesAttributes(_, _, _, maxPollingInterval)
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
