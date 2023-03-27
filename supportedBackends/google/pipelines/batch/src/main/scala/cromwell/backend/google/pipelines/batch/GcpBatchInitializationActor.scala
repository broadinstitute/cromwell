package cromwell.backend.google.pipelines.batch

import akka.actor.ActorRef
import com.google.api.services.lifesciences.v2beta.CloudLifeSciencesScopes
import com.google.api.services.cloudkms.v1.{CloudKMS, CloudKMSScopes}
import com.google.api.services.cloudkms.v1.model.EncryptRequest
import cromwell.core.WorkflowOptions
import scala.concurrent.Future
import com.google.auth.Credentials
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.backend.standard.{StandardInitializationActor, StandardInitializationActorParams, StandardValidatedRuntimeAttributesBuilder}
import cromwell.backend.google.pipelines.batch.GcpBatchInitializationActor._
import com.google.auth.oauth2.OAuth2Credentials
import com.google.auth.http.HttpCredentialsAdapter
import cromwell.core.io.AsyncIoActorClient
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import wom.graph.CommandCallNode
import com.google.api.services.storage.StorageScopes
import cromwell.filesystems.gcs.GoogleUtil._
import com.google.api.services.genomics.v2alpha1.GenomicsScopes

import cromwell.cloudsupport.gcp.auth.GoogleAuthMode.{httpTransport, jsonFactory}
import cromwell.cloudsupport.gcp.auth.{GoogleAuthMode, UserServiceAccountMode}
import org.apache.commons.codec.binary.Base64
import spray.json.{JsObject, JsString}

case class GcpBatchInitializationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  ioActor: ActorRef,
  calls: Set[CommandCallNode],
  batchConfiguration: GcpBatchConfiguration,
  serviceRegistryActor: ActorRef,
  restarting: Boolean
) extends StandardInitializationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = batchConfiguration.configurationDescriptor

}
class GcpBatchInitializationActor(batchParams: GcpBatchInitializationActorParams) extends StandardInitializationActor(batchParams) with AsyncIoActorClient {

  override lazy val ioActor: ActorRef = batchParams.ioActor
  protected val gcpBatchConfiguration: GcpBatchConfiguration = batchParams.batchConfiguration
  protected val workflowOptions: WorkflowOptions = workflowDescriptor.workflowOptions


  //private lazy val ioEc = context.system.dispatchers.lookup(Dispatcher.IoDispatcher)

  // Credentials object for the GCS API
  private lazy val gcsCredentials: Future[Credentials] = gcpBatchConfiguration
    .batchAttributes
    .auths
    .gcs
    .retryCredentials(workflowOptions, List(StorageScopes
      .DEVSTORAGE_FULL_CONTROL))

  // Credentials object for the Genomics API
  private lazy val genomicsCredentials: Future[Credentials] = gcpBatchConfiguration.batchAttributes.auths.genomics
                                                                                    .retryCredentials(workflowOptions, List(
                                                                                      CloudLifeSciencesScopes
                                                                                        .CLOUD_PLATFORM,
                                                                                      GenomicsScopes.GENOMICS
                                                                                    ))

  val privateDockerEncryptionKeyName: Option[String] = {
    val optionsEncryptionKey = workflowOptions.get(GoogleAuthMode.DockerCredentialsEncryptionKeyNameKey).toOption
    optionsEncryptionKey.orElse(gcpBatchConfiguration.dockerEncryptionKeyName)
  }

  val privateDockerToken: Option[String] = {
    val optionsDockerToken = workflowOptions.get(GoogleAuthMode.DockerCredentialsTokenKey).toOption
    optionsDockerToken.orElse(gcpBatchConfiguration.dockerToken)
  }

  lazy val privateDockerEncryptedToken: Option[String] = {
    val effectiveAuth: Option[GoogleAuthMode] = {
      // Use the user service account if a user service account value is provided in the workflow options and there's
      // a user service account auth in the list of auths.
      // TODO is it okay that this would silently ignore user auths if there isn't one defined in the config list of auths?
      // That doesn't seem great but it's effectively what the existing code around user service accounts appears to be doing.
      val userServiceAccountAuth: Option[GoogleAuthMode] = for {
        _ <- workflowOptions.get(GoogleAuthMode.UserServiceAccountKey).toOption
        usaAuth <- gcpBatchConfiguration.googleConfig.authsByName.values collectFirst { case u: UserServiceAccountMode => u }
      } yield usaAuth

      def encryptionAuthFromConfig: Option[GoogleAuthMode] = gcpBatchConfiguration.dockerEncryptionAuthName.flatMap { name =>
        gcpBatchConfiguration.googleConfig.auth(name).toOption
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
      credentials = auth.credentials(workflowOptions.get(_).get, List(CloudKMSScopes.CLOUD_PLATFORM))
      encrypted = encryptKms(key, credentials, plain)
    } yield encrypted
  }

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    GcpBatchRuntimeAttributes
      .runtimeAttributesBuilder(gcpBatchConfiguration)

  override lazy val workflowPaths: Future[GcpBatchWorkflowPaths] = for {
    gcsCred <- gcsCredentials
    genomicsCred <- genomicsCredentials
    validatedPathBuilders <- pathBuilders
  } yield new GcpBatchWorkflowPaths(
      workflowDescriptor, gcsCred, genomicsCred, gcpBatchConfiguration, validatedPathBuilders)


  override lazy val initializationData: Future[GcpBackendInitializationData] = for {
    batchWorkflowPaths <- workflowPaths
    //vpcNetworkAndSubnetworkProjectLabels <- vpcNetworkAndSubnetworkProjectLabelsFuture()
  } yield GcpBackendInitializationData(
    workflowPaths = batchWorkflowPaths,
    runtimeAttributesBuilder = runtimeAttributesBuilder,
    gcpBatchConfiguration = gcpBatchConfiguration,
    privateDockerEncryptionKeyName = privateDockerEncryptionKeyName,
    privateDockerEncryptedToken = privateDockerEncryptedToken
  )
  //add in gcs credentials if necessary

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    for {
      paths <- workflowPaths
      _ = publishWorkflowRoot(paths.workflowRoot.pathAsString)
      data <- initializationData
    } yield Option(data)
  }

  override lazy val ioCommandBuilder: GcsBatchCommandBuilder.type = GcsBatchCommandBuilder

}

object GcpBatchInitializationActor {

  def encryptKms(keyName: String, credentials: OAuth2Credentials, plainText: String): String = {
    val httpCredentialsAdapter = new HttpCredentialsAdapter(credentials)
    val kms = new CloudKMS.Builder(httpTransport, jsonFactory, httpCredentialsAdapter)
      .setApplicationName("cromwell")
      .build()

    val request = new EncryptRequest().encodePlaintext(plainText.toCharArray.map(_.toByte))
    val response = kms.projects.locations.keyRings.cryptoKeys.encrypt(keyName, request).execute
    response.getCiphertext
  }

}
