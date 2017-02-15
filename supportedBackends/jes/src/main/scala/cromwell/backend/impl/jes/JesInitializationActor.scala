package cromwell.backend.impl.jes

import java.io.IOException

import akka.actor.ActorRef
import com.google.cloud.storage.contrib.nio.CloudStorageOptions
import cromwell.backend.impl.jes.authentication.{GcsLocalizing, JesAuthInformation}
import cromwell.backend.standard.{StandardInitializationActor, StandardInitializationActorParams, StandardValidatedRuntimeAttributesBuilder}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.core.io.AsyncIo
import cromwell.filesystems.gcs.GcsBatchCommandBuilder
import cromwell.filesystems.gcs.auth.{ClientSecrets, GoogleAuthMode}
import spray.json.JsObject
import wdl4s.TaskCall

import scala.concurrent.Future

case class JesInitializationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  ioActor: ActorRef,
  calls: Set[TaskCall],
  jesConfiguration: JesConfiguration,
  serviceRegistryActor: ActorRef
) extends StandardInitializationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = jesConfiguration.configurationDescriptor
}

class JesInitializationActor(jesParams: JesInitializationActorParams)
  extends StandardInitializationActor(jesParams) with AsyncIo with GcsBatchCommandBuilder {

  override lazy val ioActor = jesParams.ioActor
  private val jesConfiguration = jesParams.jesConfiguration
  
  context.become(ioReceive orElse receive)

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    JesRuntimeAttributes.runtimeAttributesBuilder(jesConfiguration)

  private[jes] lazy val refreshTokenAuth: Option[JesAuthInformation] = {
    for {
      clientSecrets <- List(jesConfiguration.jesAttributes.auths.gcs) collectFirst { case s: ClientSecrets => s }
      token <- workflowDescriptor.workflowOptions.get(GoogleAuthMode.RefreshTokenOptionKey).toOption
    } yield GcsLocalizing(clientSecrets, token)
  }

  private lazy val genomics = jesConfiguration.genomicsFactory.withOptions(workflowDescriptor.workflowOptions)
  // FIXME: workflow paths indirectly re create part of those credentials via the GcsPathBuilder
  // This is unnecessary duplication of credentials. They are needed here so they can be added to the initialization data
  // and used to retrieve docker hashes
  private lazy val gcsCredentials = jesConfiguration.jesAuths.gcs.credential(workflowDescriptor.workflowOptions)

  override lazy val workflowPaths: JesWorkflowPaths =
    new JesWorkflowPaths(workflowDescriptor, jesConfiguration)(context.system)

  override lazy val initializationData: JesBackendInitializationData =
        JesBackendInitializationData(workflowPaths, runtimeAttributesBuilder, jesConfiguration, gcsCredentials, genomics)

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    publishWorkflowRoot(workflowPaths.workflowRoot.pathAsString)
    if (jesConfiguration.needAuthFileUpload) {
      writeAuthenticationFile(workflowPaths) map { _ => Option(initializationData) } recoverWith {
        case failure => Future.failed(new IOException("Failed to upload authentication file", failure)) 
      }
    } else {
      Future.successful(Option(initializationData))
    }
  }

  private def writeAuthenticationFile(workflowPath: JesWorkflowPaths): Future[Unit] = {
    generateAuthJson(jesConfiguration.dockerCredentials, refreshTokenAuth) map { content =>
      val path = workflowPath.gcsAuthFilePath
      workflowLogger.info(s"Creating authentication file for workflow ${workflowDescriptor.id} at \n $path")
      writeAsync(path, content, Seq(CloudStorageOptions.withMimeType("application/json")))
    } getOrElse Future.successful(())
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
