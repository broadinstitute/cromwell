package cromwell.backend.impl.jes

import java.io.IOException
import java.nio.file.StandardOpenOption

import akka.actor.ActorRef
import akka.http.scaladsl.model.StatusCodes
import com.google.api.services.genomics.Genomics
import com.google.auth.Credentials
import com.google.cloud.storage.StorageException
import com.google.cloud.storage.contrib.nio.CloudStorageOptions
import cromwell.backend.impl.jes.authentication.{GcsLocalizing, JesAuthObject, JesDockerCredentials}
import cromwell.backend.standard.{StandardInitializationActor, StandardInitializationActorParams, StandardValidatedRuntimeAttributesBuilder}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.core.CromwellFatalException
import cromwell.core.io.AsyncIo
import cromwell.filesystems.gcs.auth.{ClientSecrets, GoogleAuthMode}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import spray.json.{JsObject, JsTrue}
import wdl4s.wdl.WdlTaskCall

import scala.concurrent.Future

case class JesInitializationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  ioActor: ActorRef,
  calls: Set[WdlTaskCall],
  jesConfiguration: JesConfiguration,
  serviceRegistryActor: ActorRef,
  restarting: Boolean
) extends StandardInitializationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = jesConfiguration.configurationDescriptor
}

class JesInitializationActor(jesParams: JesInitializationActorParams)
  extends StandardInitializationActor(jesParams) with AsyncIo with GcsBatchCommandBuilder {

  override lazy val ioActor = jesParams.ioActor
  private val jesConfiguration = jesParams.jesConfiguration
  private val workflowOptions = workflowDescriptor.workflowOptions
  implicit private val system = context.system

  context.become(ioReceive orElse receive)

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    JesRuntimeAttributes.runtimeAttributesBuilder(jesConfiguration)

  // From the gcs auth and the workflow options, optionally builds a GcsLocalizing that contains
  // the information (client Id/Secrets + refresh token) that will be uploaded to Gcs before the workflow start
  private[jes] lazy val refreshTokenAuth: Option[JesAuthObject] = {
    for {
      clientSecrets <- List(jesConfiguration.jesAttributes.auths.gcs) collectFirst { case s: ClientSecrets => s }
      token <- workflowDescriptor.workflowOptions.get(GoogleAuthMode.RefreshTokenOptionKey).toOption
    } yield GcsLocalizing(clientSecrets, token)
  }

  // Credentials object for the GCS API
  private lazy val gcsCredentials: Future[Credentials] =
    jesConfiguration.jesAttributes.auths.gcs.credential(workflowOptions)

  // Credentials object for the Genomics API
  private lazy val genomicsCredentials: Future[Credentials] =
    jesConfiguration.jesAttributes.auths.genomics.credential(workflowOptions)

  // Genomics object to access the Genomics API
  private lazy val genomics: Future[Genomics] = {
    genomicsCredentials map jesConfiguration.genomicsFactory.fromCredentials
  }

  override lazy val workflowPaths: Future[JesWorkflowPaths] = for {
    gcsCred <- gcsCredentials
    genomicsCred <- genomicsCredentials
  } yield new JesWorkflowPaths(workflowDescriptor, gcsCred, genomicsCred, jesConfiguration)

  override lazy val initializationData: Future[JesBackendInitializationData] = for {
    jesWorkflowPaths <- workflowPaths
    gcsCreds <- gcsCredentials
    genomicsFactory <- genomics
  } yield JesBackendInitializationData(jesWorkflowPaths, runtimeAttributesBuilder, jesConfiguration, gcsCreds, genomicsFactory)

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    def fileUpload(paths: JesWorkflowPaths): Future[Unit] = {
      writeAuthenticationFile(paths, jesConfiguration.jesAttributes.restrictMetadataAccess, jesConfiguration.dockerCredentials) recoverWith {
        case CromwellFatalException(e: StorageException) if e.getCode == StatusCodes.PreconditionFailed.intValue && jesParams.restarting =>
          workflowLogger.debug(s"Authentication file already exists but this is a restart, proceeding.")
          Future.successful(())
        case CromwellFatalException(e: StorageException) if e.getCode == StatusCodes.PreconditionFailed.intValue =>
            Future.failed(new IOException(s"Failed to upload authentication file:" +
              " there was already a file at the same location and this workflow was not being restarted."))
        case failure => Future.failed(new IOException(s"Failed to upload authentication file", failure))
      }
    }

    for {
      paths <- workflowPaths
      _ = publishWorkflowRoot(paths.workflowRoot.pathAsString)
      _ <- if (jesConfiguration.needAuthFileUpload) fileUpload(paths) else Future.successful(())
      data <- initializationData
    } yield Option(data)
  }

  private def writeAuthenticationFile(workflowPath: JesWorkflowPaths,
                                      restrictMetadataAccess: Boolean,
                                      dockerCredentials: Option[JesDockerCredentials]): Future[Unit] = {
    val authObjects = List(dockerCredentials, refreshTokenAuth).flatten
    generateAuthJson(authObjects, restrictMetadataAccess) map { content =>
      val path = workflowPath.gcsAuthFilePath
      workflowLogger.info(s"Creating authentication file for workflow ${workflowDescriptor.id} at \n $path")
      val openOptions = Seq(
        CloudStorageOptions.withMimeType("application/json"),
        // Will fail if the file already exists - 
        // In case of a restart this avoids rewriting the file as the JES job might be trying to localize it at the same time 
        StandardOpenOption.CREATE_NEW
      )
      writeAsync(path, content, openOptions)
    } getOrElse Future.successful(())
  }

  def generateAuthJson(authObjects: List[JesAuthObject], restrictMetadataAccess: Boolean): Option[String] = {
    def generateAuthObject(): Map[String, JsObject] = {
      if (authObjects.nonEmpty) {
        val authObjectMaps = authObjects map { _.toMap }
        Map("auths" -> JsObject(authObjectMaps.reduce(_ ++ _) map { case (k, v) => k -> JsObject.apply(v) }))
      } else Map.empty
    }

    val authMap = generateAuthObject()
    val jsonMap = if (restrictMetadataAccess) authMap ++ Map("restrictMetadataAccess" -> JsTrue) else authMap

    if (jsonMap.nonEmpty) Option(JsObject(jsonMap).prettyPrint)
    else None
  }
}
