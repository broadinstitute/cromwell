package cromwell.backend.google.pipelines.common

import java.io.IOException

import akka.actor.ActorRef
import com.google.auth.Credentials
import com.google.cloud.storage.contrib.nio.CloudStorageOptions
import cromwell.backend.google.pipelines.common.JesInitializationActor.AuthFileAlreadyExistsException
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory
import cromwell.backend.google.pipelines.common.authentication.{GcsLocalizing, JesAuthObject, JesDockerCredentials}
import cromwell.backend.standard.{StandardInitializationActor, StandardInitializationActorParams, StandardValidatedRuntimeAttributesBuilder}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.cloudsupport.gcp.auth.{ClientSecrets, GoogleAuthMode}
import cromwell.core.io.AsyncIoActorClient
import cromwell.core.path.Path
import cromwell.filesystems.gcs.GoogleUtil._
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import spray.json.{JsObject, JsTrue}
import wom.graph.CommandCallNode

import scala.concurrent.Future

case class JesInitializationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  ioActor: ActorRef,
  calls: Set[CommandCallNode],
  jesConfiguration: PipelinesApiConfiguration,
  serviceRegistryActor: ActorRef,
  restarting: Boolean
) extends StandardInitializationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = jesConfiguration.configurationDescriptor
}

object JesInitializationActor {
  private case class AuthFileAlreadyExistsException(path: Path) extends IOException(s"Failed to upload authentication file at $path:" +
    s" there was already a file at the same location and this workflow was not being restarted.")
}

class JesInitializationActor(jesParams: JesInitializationActorParams)
  extends StandardInitializationActor(jesParams) with AsyncIoActorClient {

  override lazy val ioActor = jesParams.ioActor
  private val jesConfiguration = jesParams.jesConfiguration
  private val workflowOptions = workflowDescriptor.workflowOptions

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    PipelinesApiRuntimeAttributes.runtimeAttributesBuilder(jesConfiguration)

  // From the gcs auth and the workflow options, optionally builds a GcsLocalizing that contains
  // the information (client Id/Secrets + refresh token) that will be uploaded to Gcs before the workflow start
  private[pipelines] lazy val refreshTokenAuth: Option[JesAuthObject] = {
    for {
      clientSecrets <- List(jesConfiguration.jesAttributes.auths.gcs) collectFirst { case s: ClientSecrets => s }
      token <- workflowDescriptor.workflowOptions.get(GoogleAuthMode.RefreshTokenOptionKey).toOption
    } yield GcsLocalizing(clientSecrets, token)
  }

  // Credentials object for the GCS API
  private lazy val gcsCredentials: Future[Credentials] =
    jesConfiguration.jesAttributes.auths.gcs.retryCredential(workflowOptions)

  // Credentials object for the Genomics API
  private lazy val genomicsCredentials: Future[Credentials] =
    jesConfiguration.jesAttributes.auths.genomics.retryCredential(workflowOptions)

  // Genomics object to access the Genomics API
  private lazy val genomics: Future[PipelinesApiRequestFactory] = {
    genomicsCredentials map jesConfiguration.genomicsFactory.fromCredentials
  }

  override lazy val workflowPaths: Future[PipelinesApiWorkflowPaths] = for {
    gcsCred <- gcsCredentials
    genomicsCred <- genomicsCredentials
    validatedPathBuilders <- pathBuilders
  } yield new PipelinesApiWorkflowPaths(workflowDescriptor, gcsCred, genomicsCred, jesConfiguration, validatedPathBuilders)

  override lazy val initializationData: Future[PipelinesApiBackendInitializationData] = for {
    jesWorkflowPaths <- workflowPaths
    gcsCreds <- gcsCredentials
    genomicsFactory <- genomics
  } yield PipelinesApiBackendInitializationData(jesWorkflowPaths, runtimeAttributesBuilder, jesConfiguration, gcsCreds, genomicsFactory)

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    def fileUpload(paths: PipelinesApiWorkflowPaths) = asyncIo.existsAsync(paths.gcsAuthFilePath) flatMap {
      // if we aren't restarting then something is definitely wrong
      case true if !jesParams.restarting =>
        Future.failed(AuthFileAlreadyExistsException(paths.gcsAuthFilePath))
      case true => 
        workflowLogger.debug(s"Authentication file already exists but this is a restart, proceeding.")
        Future.successful(())
      case false =>
        writeAuthenticationFile(paths, jesConfiguration.jesAttributes.restrictMetadataAccess, jesConfiguration.dockerCredentials) recoverWith {
          case failure => Future.failed(new IOException(s"Failed to upload authentication file", failure))
        }
    } recoverWith {
      case failure => Future.failed(new IOException(s"Failed to upload authentication file", failure))
    }
    
    for {
      paths <- workflowPaths
      _ = publishWorkflowRoot(paths.workflowRoot.pathAsString)
      _ <- if (jesConfiguration.needAuthFileUpload) fileUpload(paths) else Future.successful(())
      data <- initializationData
    } yield Option(data)
  }

  private def writeAuthenticationFile(workflowPath: PipelinesApiWorkflowPaths,
                                      restrictMetadataAccess: Boolean,
                                      dockerCredentials: Option[JesDockerCredentials]): Future[Unit] = {
    val authObjects = List(dockerCredentials, refreshTokenAuth).flatten
    generateAuthJson(authObjects, restrictMetadataAccess) map { content =>
      val path = workflowPath.gcsAuthFilePath
      workflowLogger.info(s"Creating authentication file for workflow ${workflowDescriptor.id} at \n $path")
      val openOptions = List(CloudStorageOptions.withMimeType("application/json"))
      asyncIo.writeAsync(path, content, openOptions)
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

  override lazy val ioCommandBuilder = GcsBatchCommandBuilder
}
