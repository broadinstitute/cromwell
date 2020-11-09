package cromwell.backend.google.pipelines.common

import akka.actor.ActorRef
import com.google.api.client.http.{HttpRequest, HttpResponse}
import com.google.api.services.cloudkms.v1.model.EncryptRequest
import com.google.api.services.cloudkms.v1.{CloudKMS, CloudKMSScopes}
import com.google.api.services.cloudresourcemanager.CloudResourceManager
import com.google.api.services.genomics.GenomicsScopes
import com.google.api.services.lifesciences.v2beta.CloudLifeSciencesScopes
import com.google.api.services.storage.StorageScopes
import com.google.auth.Credentials
import com.google.auth.http.HttpCredentialsAdapter
import com.google.auth.oauth2.OAuth2Credentials
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.VirtualPrivateCloudConfiguration
import cromwell.backend.google.pipelines.common.PipelinesApiInitializationActor._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory
import cromwell.backend.standard.{StandardInitializationActor, StandardInitializationActorParams, StandardValidatedRuntimeAttributesBuilder}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode.{httpTransport, jsonFactory}
import cromwell.cloudsupport.gcp.auth.{GoogleAuthMode, UserServiceAccountMode}
import cromwell.core.{Dispatcher, WorkflowOptions}
import cromwell.core.io.AsyncIoActorClient
import cromwell.filesystems.gcs.GoogleUtil._
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import _root_.io.circe.Decoder
import _root_.io.circe.generic.semiauto.deriveDecoder
import _root_.io.circe.parser.decode
import org.apache.commons.codec.binary.Base64
import org.apache.commons.lang3.exception.ExceptionUtils
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

  override lazy val ioActor: ActorRef = pipelinesParams.ioActor
  protected val pipelinesConfiguration: PipelinesApiConfiguration = pipelinesParams.jesConfiguration
  protected val workflowOptions: WorkflowOptions = workflowDescriptor.workflowOptions
  private lazy val ioEc = context.system.dispatchers.lookup(Dispatcher.IoDispatcher)

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    PipelinesApiRuntimeAttributes.runtimeAttributesBuilder(pipelinesConfiguration)

  // Credentials object for the GCS API
  private lazy val gcsCredentials: Future[Credentials] = pipelinesConfiguration.papiAttributes.auths.gcs
      .retryCredentials(workflowOptions, List(StorageScopes.DEVSTORAGE_FULL_CONTROL))

  // Credentials object for the Genomics API
  private lazy val genomicsCredentials: Future[Credentials] = pipelinesConfiguration.papiAttributes.auths.genomics
    .retryCredentials(workflowOptions, List(
      CloudLifeSciencesScopes.CLOUD_PLATFORM,
      GenomicsScopes.GENOMICS,
      /*
      Genomics Pipelines API v1alpha2 requires the COMPUTE scope. Does not seem to be required for either v2alpha1 or v2beta.
       */
      GenomicsScopes.COMPUTE,
      /*
      Used to write so-called "auth" files. The `gcsAuthFilePath` could probably be refactored such that *this* created
      `genomicsCredentials` doesn't actually need DEVSTORAGE_FULL_CONTROL, but it's also not clear how the magic
      Genomics Pipelines API parameter "__extra_config_gcs_path" works nor where it's documented.

      See also:
      - cromwell.backend.google.pipelines.common.PipelinesApiWorkflowPaths#gcsAuthFilePath
      - https://github.com/broadinstitute/cromwell/pull/2435
      - cromwell.backend.google.pipelines.common.PipelinesApiConfiguration#needAuthFileUpload
      - cromwell.cloudsupport.gcp.auth.GoogleAuthMode#requiresAuthFile
      - cromwell.backend.google.pipelines.common.PipelinesApiAsyncBackendJobExecutionActor#gcsAuthParameter
       */
      StorageScopes.DEVSTORAGE_FULL_CONTROL,
    ))

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
      credentials = auth.credentials(workflowOptions.get(_).get, List(CloudKMSScopes.CLOUD_PLATFORM))
      encrypted = encryptKms(key, credentials, plain)
    } yield encrypted
  }

  private def vpcNetworkAndSubnetworkProjectLabelsFuture(): Future[Option[VpcAndSubnetworkProjectLabelValues]] = {
    def googleProject(descriptor: BackendWorkflowDescriptor): String = {
      descriptor.workflowOptions.getOrElse(WorkflowOptionKeys.GoogleProject, pipelinesParams.jesConfiguration.papiAttributes.project)
    }

    def projectMetadataRequest(vpcConfig: VirtualPrivateCloudConfiguration): Future[HttpRequest] = {
      Future {
        val credentials = vpcConfig.auth.credentials(workflowOptions.get(_).get, List(CloudLifeSciencesScopes.CLOUD_PLATFORM))

        val httpCredentialsAdapter = new HttpCredentialsAdapter(credentials)
        val cloudResourceManagerBuilder = new CloudResourceManager
        .Builder(GoogleAuthMode.httpTransport, GoogleAuthMode.jsonFactory, httpCredentialsAdapter)
          .setApplicationName(pipelinesConfiguration.googleConfig.applicationName)
          .build()

        val project = cloudResourceManagerBuilder.projects().get(googleProject(workflowDescriptor))

        project.buildHttpRequest()
      }
    }

    def projectMetadataResponseToLabels(httpResponse: HttpResponse): Future[ProjectLabels] = {
      implicit val googleProjectMetadataLabelDecoder: Decoder[ProjectLabels] = deriveDecoder
      Future.fromTry(decode[ProjectLabels](httpResponse.parseAsString()).toTry).recoverWith {
        case e: Throwable => Future.failed(new RuntimeException(s"Failed to parse labels from project metadata response from Google Cloud Resource Manager API. " +
          s"${ExceptionUtils.getMessage(e)}", e))
      }
    }

    def networkLabelsFromProjectLabels(vpcConfig: VirtualPrivateCloudConfiguration, projectLabels: ProjectLabels): Option[VpcAndSubnetworkProjectLabelValues] = {
      projectLabels.labels.get(vpcConfig.name) map { vpcNetworkLabelValue =>
        val subnetworkLabelOption = vpcConfig.subnetwork.flatMap { s =>
          projectLabels.labels.collectFirst {
            case (labelName, labelValue) if labelName.equals(s) => labelValue
          }
        }

        VpcAndSubnetworkProjectLabelValues(vpcNetworkLabelValue, subnetworkLabelOption)
      }
    }

    def fetchVpcLabelsFromProjectMetadata(vpcConfig: VirtualPrivateCloudConfiguration): Future[Option[VpcAndSubnetworkProjectLabelValues]] = {
      for {
        projectMetadataResponse <- projectMetadataRequest(vpcConfig).map(_.executeAsync().get())
        projectLabels <- projectMetadataResponseToLabels(projectMetadataResponse)
      } yield networkLabelsFromProjectLabels(vpcConfig, projectLabels)
    }

    pipelinesConfiguration.papiAttributes.virtualPrivateCloudConfiguration match {
      case None => Future.successful(None)
      case Some(vpcConfig) => fetchVpcLabelsFromProjectMetadata(vpcConfig)
    }
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
    vpcNetworkAndSubnetworkProjectLabels <- vpcNetworkAndSubnetworkProjectLabelsFuture()
  } yield PipelinesApiBackendInitializationData(
    workflowPaths = jesWorkflowPaths,
    runtimeAttributesBuilder = runtimeAttributesBuilder,
    papiConfiguration = pipelinesConfiguration,
    gcsCredentials = gcsCreds,
    genomicsRequestFactory = genomicsFactory,
    privateDockerEncryptionKeyName = privateDockerEncryptionKeyName,
    privateDockerEncryptedToken = privateDockerEncryptedToken,
    vpcNetworkAndSubnetworkProjectLabels = vpcNetworkAndSubnetworkProjectLabels)

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

  override lazy val ioCommandBuilder: GcsBatchCommandBuilder.type = GcsBatchCommandBuilder
}

object PipelinesApiInitializationActor {
  // For metadata publishing purposes default to using the name of a standard stream as the stream's filename.
  def defaultStandardStreamNameToFileNameMetadataMapper(pipelinesApiJobPaths: PipelinesApiJobPaths, streamName: String): String = streamName

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
