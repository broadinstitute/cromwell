package cromwell.backend.google.pipelines.v1alpha2

import java.io.IOException

import com.google.cloud.storage.contrib.nio.CloudStorageOptions
import cromwell.backend.google.pipelines.common.authentication.{GcsLocalizing, PipelinesApiAuthObject, PipelinesApiDockerCredentials}
import cromwell.backend.google.pipelines.common.{PipelinesApiInitializationActorParams, PipelinesApiJobPaths, PipelinesApiWorkflowPaths}
import cromwell.backend.google.pipelines.v1alpha2.PipelinesApiInitializationActor.AuthFileAlreadyExistsException
import cromwell.cloudsupport.gcp.auth.{ClientSecrets, GoogleAuthMode}
import cromwell.core.path.Path
import spray.json.{JsObject, JsTrue}

import scala.concurrent.Future

class PipelinesApiInitializationActor(pipelinesParams: PipelinesApiInitializationActorParams)
  extends cromwell.backend.google.pipelines.common.PipelinesApiInitializationActor(pipelinesParams) {

  // From the gcs auth and the workflow options, optionally builds a GcsLocalizing that contains
  // the information (client Id/Secrets + refresh token) that will be uploaded to Gcs before the workflow start
  private[pipelines] lazy val refreshTokenAuth: Option[PipelinesApiAuthObject] = {
    for {
      clientSecrets <- List(pipelinesConfiguration.papiAttributes.auths.gcs) collectFirst { case s: ClientSecrets => s }
      token <- workflowDescriptor.workflowOptions.get(GoogleAuthMode.RefreshTokenOptionKey).toOption
    } yield GcsLocalizing(clientSecrets, token)
  }

  // V1 needs to upload the auth file at the beginning of the workflow, so override the default beforeAll
  override def beforeAll() = {
    def fileUpload(paths: PipelinesApiWorkflowPaths) = asyncIo.existsAsync(paths.gcsAuthFilePath) flatMap {
      // if we aren't restarting then something is definitely wrong
      case true if !pipelinesParams.restarting =>
        Future.failed(AuthFileAlreadyExistsException(paths.gcsAuthFilePath))
      case true =>
        workflowLogger.debug(s"Authentication file already exists but this is a restart, proceeding.")
        Future.successful(())
      case false =>
        writeAuthenticationFile(paths, pipelinesConfiguration.papiAttributes.restrictMetadataAccess, pipelinesConfiguration.dockerCredentials) recoverWith {
          case failure => Future.failed(new IOException(s"Failed to upload authentication file", failure))
        }
    } recoverWith {
      case failure => Future.failed(new IOException(s"Failed to upload authentication file", failure))
    }

    for {
      paths <- workflowPaths
      _ <- if (pipelinesConfiguration.needAuthFileUpload) fileUpload(paths) else Future.successful(())
      data <- super.beforeAll()
    } yield data
  }

  private def writeAuthenticationFile(workflowPath: PipelinesApiWorkflowPaths,
                                      restrictMetadataAccess: Boolean,
                                      dockerCredentials: Option[PipelinesApiDockerCredentials]): Future[Unit] = {
    val authObjects = List(dockerCredentials, refreshTokenAuth).flatten
    generateAuthJson(authObjects, restrictMetadataAccess) map { content =>
      val path = workflowPath.gcsAuthFilePath
      workflowLogger.info(s"Creating authentication file for workflow ${workflowDescriptor.id} at \n $path")
      val openOptions = List(CloudStorageOptions.withMimeType("application/json"))
      asyncIo.writeAsync(path, content, openOptions)
    } getOrElse Future.successful(())
  }

  def generateAuthJson(authObjects: List[PipelinesApiAuthObject], restrictMetadataAccess: Boolean): Option[String] = {
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

  /** Overridden for v1 to account for the way the PAPI v1 controller names these files. PAPI v1 completely manages the periodic
    * delocalization of these files including their naming. */
  override def standardStreamNameToFileNameMetadataMapper(pipelinesApiJobPaths: PipelinesApiJobPaths, streamName: String): String =
    s"${pipelinesApiJobPaths.jesLogBasename}-$streamName.log"
}

object PipelinesApiInitializationActor {
  private case class AuthFileAlreadyExistsException(path: Path) extends IOException(s"Failed to upload authentication file at $path:" +
    s" there was already a file at the same location and this workflow was not being restarted.")
}
