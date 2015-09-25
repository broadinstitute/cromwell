package cromwell.engine.backend.jes

import java.net.URL

import com.typesafe.config.ConfigFactory
import cromwell.util.ConfigUtil._
import cromwell.util.ReferenceConfiguration

import scala.collection.JavaConversions._
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

  private val requiredConfig = ConfigFactory.parseMap(Map(
    "applicationName" -> "",
    "project" -> "",
    "baseExecutionBucket" -> "",
    "endpointUrl" -> "",
    "authenticationMode" -> ""
  ))

  private val optionalConfig = Option(ConfigFactory.parseMap(Map(
    "dockerAccount" -> "",
    "dockerToken" -> ""
  )))

  val referenceJesConf = ReferenceConfiguration(requiredConfig, optionalConfig, "Jes")

  def apply(): JesAttributes = {
    val jesConf = ConfigFactory.load.getConfig("backend").getConfig("jes")

    jesConf.checkValidWithWarnings(referenceJesConf)

    val applicationName = jesConf.getString("applicationName")
    val project = jesConf.getString("project")
    val executionBucket = jesConf.getString("baseExecutionBucket")

    // values requiring extended validation
    val endpointUrl: ValidationNel[String, URL] = jesConf.validateURL("endpointUrl")
    val authMode: ValidationNel[String, GcsAuthMode] = {
      jesConf.getString("authenticationMode").validateAny[GcsAuthMode, IllegalArgumentException] { GcsAuthMode.fromString }
    }

    val dockerCredentials = for {
      account <- jesConf.getStringOption("dockerAccount")
      token <- jesConf.getStringOption("dockerToken")
    } yield DockerCredentials(account, token)

    (endpointUrl |@| authMode) { JesAttributes(applicationName, project, executionBucket, _, _, dockerCredentials) } match {
      case Success(r) => r
      case Failure(f) =>
        val errorMessages = f.list.mkString(", ")
        throw new IllegalArgumentException(s"Jes Configuration is not valid: Errors: $errorMessages")
    }
  }

}
