package cromwell.engine.backend.jes

import java.net.URL

import com.typesafe.config.ConfigFactory
import cromwell.util.ConfigUtil._

import scala.language.postfixOps
import scalaz.Scalaz._
import scalaz._

case class JesAttributes(applicationName: String,
                         project: String,
                         executionBucket: String,
                         endpointUrl: URL,
                         authMode: GcsAuthMode,
                         dockerCredentials: Option[DockerCredentials])
object JesAttributes {

  private val keys = Set(
    "applicationName",
    "project",
    "baseExecutionBucket",
    "endpointUrl",
    "authenticationMode",
    "maximumPollingInterval",
    "dockerAccount",
    "dockerToken"
  )

  private val context = "Jes"

  def apply(): JesAttributes = {
    val jesConf = ConfigFactory.load.getConfig("backend").getConfig("jes")

    jesConf.warnNotRecognized(keys, context)

    val applicationName: ValidationNel[String, String] = jesConf.validateString("applicationName")
    val project: ValidationNel[String, String] = jesConf.validateString("project")
    val executionBucket: ValidationNel[String, String] = jesConf.validateString("baseExecutionBucket")
    val endpointUrl: ValidationNel[String, URL] = jesConf.validateURL("endpointUrl")
    val authMode: ValidationNel[String, GcsAuthMode] = jesConf.getString("authenticationMode") validateAny GcsAuthMode.fromString

    val dockerCredentials = for {
      account <- jesConf.getStringOption("dockerAccount")
      token <- jesConf.getStringOption("dockerToken")
    } yield DockerCredentials(account, token)

    (applicationName |@| project |@| executionBucket |@| endpointUrl |@| authMode) {
      JesAttributes(_, _, _, _, _, dockerCredentials)
    } match {
      case Success(r) => r
      case Failure(f) =>
        val errorMessages = f.list.mkString(", ")
        throw new IllegalArgumentException(s"Jes Configuration is not valid: Errors: $errorMessages")
    }
  }

}

trait JesConfiguration {
  def getAttributes: JesAttributes = JesAttributes()
}