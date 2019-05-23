package cromwell.backend.google.pipelines.common

import akka.actor.ActorRef
import com.google.auth.Credentials
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory
import cromwell.backend.standard.{StandardInitializationActor, StandardInitializationActorParams, StandardValidatedRuntimeAttributesBuilder}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.cloudsupport.gcp.auth.{GoogleAuthMode, UserServiceAccountMode}
import cromwell.core.Dispatcher
import cromwell.core.io.AsyncIoActorClient
import cromwell.filesystems.gcs.GoogleUtil._
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import org.apache.commons.codec.binary.Base64
import spray.json.{JsObject, JsString}
import wom.graph.CommandCallNode

import scala.concurrent.Future

case class PipelinesApiInitializationActorParams
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

class PipelinesApiInitializationActor(pipelinesParams: PipelinesApiInitializationActorParams)
  extends StandardInitializationActor(pipelinesParams) with AsyncIoActorClient {

  override lazy val ioActor = pipelinesParams.ioActor
  protected val pipelinesConfiguration = pipelinesParams.jesConfiguration
  protected val workflowOptions = workflowDescriptor.workflowOptions
  private lazy val ioEc = context.system.dispatchers.lookup(Dispatcher.IoDispatcher)

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    PipelinesApiRuntimeAttributes.runtimeAttributesBuilder(pipelinesConfiguration)

  // Credentials object for the GCS API
  private lazy val gcsCredentials: Future[Credentials] =
    pipelinesConfiguration.papiAttributes.auths.gcs.retryPipelinesApiCredentials(workflowOptions)

  // Credentials object for the Genomics API
  private lazy val genomicsCredentials: Future[Credentials] =
    pipelinesConfiguration.papiAttributes.auths.genomics.retryPipelinesApiCredentials(workflowOptions)

  // Genomics object to access the Genomics API
  private lazy val genomics: Future[PipelinesApiRequestFactory] = {
    genomicsCredentials map pipelinesConfiguration.genomicsFactory.fromCredentials
  }

  val privateDockerEncryptionKeyName: Option[String] = {
    val optionsEncryptionKey = workflowOptions.get(GoogleAuthMode.DockerCredentialsEncryptionKeyNameKey).toOption
    optionsEncryptionKey.orElse(pipelinesConfiguration.dockerEncryptionKeyName)
  }

  val privateDockerToken: Option[String] = {
    val optionsDockerToken = workflowOptions.get(GoogleAuthMode.DockerCredentialsTokenKey).toOption
    optionsDockerToken.orElse(pipelinesConfiguration.dockerToken)
  }

  lazy val privateDockerEncryptedToken: Option[String] = {
    val effectiveAuth: Option[GoogleAuthMode] = {
      // Use the user service account if a user service account value is provided in the workflow options and there's
      // a user service account auth in the list of auths.
      // TODO is it okay that this would silently ignore user auths if there isn't one defined in the config list of auths?
      // That doesn't seem great but it's effectively what the existing code around user service accounts appears to be doing.
      val userServiceAccountAuth: Option[GoogleAuthMode] = for {
        _ <- workflowOptions.get(GoogleAuthMode.UserServiceAccountKey).toOption
        usaAuth <- pipelinesConfiguration.googleConfig.authsByName.values collectFirst { case u: UserServiceAccountMode => u }
      } yield usaAuth

      def encryptionAuthFromConfig: Option[GoogleAuthMode] = pipelinesConfiguration.dockerEncryptionAuthName.flatMap { name =>
        pipelinesConfiguration.googleConfig.auth(name).toOption
      }
      // If there's no user service account auth in the workflow options fall back to an auth specified in config.
      userServiceAccountAuth orElse encryptionAuthFromConfig
    }

    val unencrypted: Option[String] = privateDockerToken flatMap { dockerToken =>
      new String(Base64.decodeBase64(dockerToken)).split(':') match {
        case Array(username, password) =>
          // unencrypted tokens are base64-encoded username:password
          Option(JsObject(
            Map(
              "username" -> JsString(username),
              "password" -> JsString(password)
            )).compactPrint)
        case _ => throw new RuntimeException(s"provided dockerhub token '$dockerToken' is not a base64-encoded username:password")
      }
    }

    for {
      plain <- unencrypted
      auth <- effectiveAuth
      key <- privateDockerEncryptionKeyName
      cred <- auth.apiClientGoogleCredential(k => workflowOptions.get(k).get)
    } yield GoogleAuthMode.encryptKms(key, cred, plain)
  }

  override lazy val workflowPaths: Future[PipelinesApiWorkflowPaths] = for {
    gcsCred <- gcsCredentials
    genomicsCred <- genomicsCredentials
    validatedPathBuilders <- pathBuilders
  } yield new PipelinesApiWorkflowPaths(
    workflowDescriptor, gcsCred, genomicsCred, pipelinesConfiguration, validatedPathBuilders, standardStreamNameToFileNameMetadataMapper)(ioEc)

  override lazy val initializationData: Future[PipelinesApiBackendInitializationData] = for {
    jesWorkflowPaths <- workflowPaths
    gcsCreds <- gcsCredentials
    genomicsFactory <- genomics
  } yield PipelinesApiBackendInitializationData(
    workflowPaths = jesWorkflowPaths,
    runtimeAttributesBuilder = runtimeAttributesBuilder,
    papiConfiguration = pipelinesConfiguration,
    gcsCredentials = gcsCreds,
    genomicsRequestFactory = genomicsFactory,
    privateDockerEncryptionKeyName = privateDockerEncryptionKeyName,
    privateDockerEncryptedToken = privateDockerEncryptedToken)

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    for {
      paths <- workflowPaths
      _ = publishWorkflowRoot(paths.workflowRoot.pathAsString)
      // Validate the google-labels workflow options, and only succeed initialization if they're good:
      _ <- Future.fromTry(GoogleLabels.fromWorkflowOptions(workflowOptions))
      data <- initializationData
    } yield Option(data)
  }

  def standardStreamNameToFileNameMetadataMapper(pipelinesApiJobPaths: PipelinesApiJobPaths, streamName: String): String =
    PipelinesApiInitializationActor.defaultStandardStreamNameToFileNameMetadataMapper(pipelinesApiJobPaths, streamName)

  override lazy val ioCommandBuilder = GcsBatchCommandBuilder
}

object PipelinesApiInitializationActor {
  // For metadata publishing purposes default to using the name of a standard stream as the stream's filename.
  def defaultStandardStreamNameToFileNameMetadataMapper(pipelinesApiJobPaths: PipelinesApiJobPaths, streamName: String): String = streamName
}
