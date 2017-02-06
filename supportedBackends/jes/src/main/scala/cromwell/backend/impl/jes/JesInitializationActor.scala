package cromwell.backend.impl.jes

import java.io.IOException

import akka.actor.ActorRef
import cromwell.backend.impl.jes.authentication.{GcsLocalizing, JesAuthInformation}
import cromwell.backend.impl.jes.io._
import cromwell.backend.standard.{StandardInitializationActor, StandardInitializationActorParams, StandardValidatedRuntimeAttributesBuilder}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.filesystems.gcs.auth.{ClientSecrets, GoogleAuthMode}
import spray.json.JsObject
import wdl4s.TaskCall

import scala.concurrent.Future
import scala.util.Try

case class JesInitializationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  calls: Set[TaskCall],
  jesConfiguration: JesConfiguration,
  serviceRegistryActor: ActorRef
) extends StandardInitializationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = jesConfiguration.configurationDescriptor
}

class JesInitializationActor(jesParams: JesInitializationActorParams)
  extends StandardInitializationActor(jesParams) {

  private val jesConfiguration = jesParams.jesConfiguration

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    JesRuntimeAttributes.runtimeAttributesBuilder(jesConfiguration)

  private[jes] lazy val refreshTokenAuth: Option[JesAuthInformation] = {
    for {
      clientSecrets <- List(jesConfiguration.jesAttributes.auths.gcs) collectFirst { case s: ClientSecrets => s }
      token <- workflowDescriptor.workflowOptions.get(GoogleAuthMode.RefreshTokenOptionKey).toOption
    } yield GcsLocalizing(clientSecrets, token)
  }

  private lazy val genomics = jesConfiguration.genomicsFactory.withOptions(workflowDescriptor.workflowOptions)

  override lazy val workflowPaths: JesWorkflowPaths =
    new JesWorkflowPaths(workflowDescriptor, jesConfiguration)(context.system)


  override lazy val initializationData: JesBackendInitializationData =
    JesBackendInitializationData(workflowPaths, runtimeAttributesBuilder, jesConfiguration, genomics)

  override def beforeAll(): Future[Option[BackendInitializationData]] = Future.fromTry(Try {
    if (jesConfiguration.needAuthFileUpload) writeAuthenticationFile(workflowPaths)
    publishWorkflowRoot(workflowPaths.workflowRoot.toString)
    Option(initializationData)
  })

  private def writeAuthenticationFile(workflowPath: JesWorkflowPaths): Unit = {
    generateAuthJson(jesConfiguration.dockerCredentials, refreshTokenAuth) foreach { content =>
      val path = workflowPath.gcsAuthFilePath
      workflowLogger.info(s"Creating authentication file for workflow ${workflowDescriptor.id} at \n ${path.toUri}")
      try {
        path.writeAsJson(content)
      } catch {
        case exception: Exception => throw new IOException("Failed to upload authentication file", exception)
      }
    }
  }

  def generateAuthJson(authInformation: Option[JesAuthInformation]*): Option[String] = {
    authInformation.flatten map { _.toMap } match {
      case Nil => None
      case jsons =>
        val authsValues = jsons.reduce(_ ++ _) mapValues JsObject.apply
        Option(JsObject("auths" -> JsObject(authsValues)).prettyPrint)
    }
  }
}
