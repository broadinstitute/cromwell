package cromwell.backend.google.batch.actors

import _root_.io.circe.Decoder
import _root_.io.circe.generic.semiauto.deriveDecoder
import _root_.io.circe.parser.decode
import akka.actor.ActorRef
import com.google.api.client.http.{HttpRequest, HttpResponse}
import com.google.api.services.cloudkms.v1.model.EncryptRequest
import com.google.api.services.cloudkms.v1.{CloudKMS, CloudKMSScopes}
import com.google.api.services.cloudresourcemanager.{CloudResourceManager, CloudResourceManagerScopes}
import com.google.api.services.genomics.v2alpha1.GenomicsScopes
import com.google.api.services.lifesciences.v2beta.CloudLifeSciencesScopes
import com.google.api.services.storage.StorageScopes
import com.google.auth.Credentials
import com.google.auth.http.HttpCredentialsAdapter
import com.google.auth.oauth2.OAuth2Credentials
import cromwell.backend.google.batch._
import cromwell.backend.google.batch.actors.GcpBatchInitializationActor._
import cromwell.backend.google.batch.api.GcpBatchRequestFactoryImpl
import cromwell.backend.google.batch.models.GcpBatchConfigurationAttributes.{
  VirtualPrivateCloudConfiguration,
  VirtualPrivateCloudLabels,
  VirtualPrivateCloudLiterals
}
import cromwell.backend.google.batch.models._
import cromwell.backend.google.batch.runnable.WorkflowOptionKeys
import cromwell.backend.standard.{
  StandardInitializationActor,
  StandardInitializationActorParams,
  StandardValidatedRuntimeAttributesBuilder
}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode.{httpTransport, jsonFactory}
import cromwell.cloudsupport.gcp.auth.{GoogleAuthMode, UserServiceAccountMode}
import cromwell.core.io.AsyncIoActorClient
import cromwell.core.{Dispatcher, WorkflowOptions}
import cromwell.filesystems.gcs.GoogleUtil._
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import org.apache.commons.codec.binary.Base64
import org.apache.commons.lang3.exception.ExceptionUtils
import spray.json.{JsObject, JsString}
import wom.graph.CommandCallNode

import scala.concurrent.Future
import scala.util.Try
import scala.util.control.NonFatal

case class GcpBatchInitializationActorParams(
  workflowDescriptor: BackendWorkflowDescriptor,
  ioActor: ActorRef,
  calls: Set[CommandCallNode],
  batchConfiguration: GcpBatchConfiguration,
  serviceRegistryActor: ActorRef,
  restarting: Boolean
) extends StandardInitializationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = batchConfiguration.configurationDescriptor

}
class GcpBatchInitializationActor(batchParams: GcpBatchInitializationActorParams)
    extends StandardInitializationActor(batchParams)
    with AsyncIoActorClient {

  override lazy val ioActor: ActorRef = batchParams.ioActor
  protected val gcpBatchConfiguration: GcpBatchConfiguration = batchParams.batchConfiguration
  protected val workflowOptions: WorkflowOptions = workflowDescriptor.workflowOptions
  private lazy val ioEc = context.system.dispatchers.lookup(Dispatcher.IoDispatcher)

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    GcpBatchRuntimeAttributes
      .runtimeAttributesBuilder(gcpBatchConfiguration)

  // Credentials object for the GCS API
  private lazy val gcsCredentials: Future[Credentials] = gcpBatchConfiguration.batchAttributes.auths.gcs
    .retryCredentials(workflowOptions, List(StorageScopes.DEVSTORAGE_FULL_CONTROL))

  // Credentials object for the Genomics API
  private lazy val genomicsCredentials: Future[Credentials] = gcpBatchConfiguration.batchAttributes.auths.genomics
    .retryCredentials(workflowOptions,
                      List(
                        CloudLifeSciencesScopes.CLOUD_PLATFORM,
                        GenomicsScopes.GENOMICS
                      )
    )

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
        usaAuth <- gcpBatchConfiguration.googleConfig.authsByName.values collectFirst {
          case u: UserServiceAccountMode => u
        }
      } yield usaAuth

      def encryptionAuthFromConfig: Option[GoogleAuthMode] = gcpBatchConfiguration.dockerEncryptionAuthName.flatMap {
        name =>
          gcpBatchConfiguration.googleConfig.auth(name).toOption
      }
      // If there's no user service account auth in the workflow options fall back to an auth specified in config.
      userServiceAccountAuth orElse encryptionAuthFromConfig
    }

    val unencrypted: Option[String] = privateDockerToken flatMap { dockerToken =>
      new String(Base64.decodeBase64(dockerToken)).split(':') match {
        case Array(username, password) =>
          // unencrypted tokens are base64-encoded username:password
          Option(
            JsObject(
              Map(
                "username" -> JsString(username),
                "password" -> JsString(password)
              )
            ).compactPrint
          )
        case _ =>
          throw new RuntimeException(
            s"provided dockerhub token '$dockerToken' is not a base64-encoded username:password"
          )
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
    def googleProject(descriptor: BackendWorkflowDescriptor): String =
      descriptor.workflowOptions.getOrElse(WorkflowOptionKeys.GoogleProject,
                                           batchParams.batchConfiguration.batchAttributes.project
      )

    def projectMetadataRequest(vpcConfig: VirtualPrivateCloudLabels): Future[HttpRequest] =
      Future {
        val credentials = vpcConfig.auth.credentials(
          workflowOptions
            .get(_)
            .getOrElse(throw new RuntimeException("Unable to find the necessary workflow option for auth credentials")),
          List(CloudResourceManagerScopes.CLOUD_PLATFORM)
        )

        val httpCredentialsAdapter = new HttpCredentialsAdapter(credentials)
        val cloudResourceManagerBuilder = new CloudResourceManager.Builder(GoogleAuthMode.httpTransport,
                                                                           GoogleAuthMode.jsonFactory,
                                                                           httpCredentialsAdapter
        )
          .setApplicationName(gcpBatchConfiguration.googleConfig.applicationName)
          .build()

        val project = cloudResourceManagerBuilder.projects().get(googleProject(workflowDescriptor))

        project.buildHttpRequest()
      }

    def projectMetadataResponseToLabels(httpResponse: HttpResponse): Future[ProjectLabels] = {
      implicit val googleProjectMetadataLabelDecoder: Decoder[ProjectLabels] = deriveDecoder
      Future.fromTry(decode[ProjectLabels](httpResponse.parseAsString()).toTry).recoverWith { case NonFatal(e) =>
        Future.failed(
          new RuntimeException(
            s"Failed to parse labels from project metadata response from Google Cloud Resource Manager API. " +
              s"${ExceptionUtils.getMessage(e)}",
            e
          )
        )
      }
    }

    def networkLabelsFromProjectLabels(vpcConfig: VirtualPrivateCloudLabels,
                                       projectLabels: ProjectLabels
    ): Option[VpcAndSubnetworkProjectLabelValues] =
      projectLabels.labels.get(vpcConfig.network) map { vpcNetworkLabelValue =>
        val subnetworkLabelOption = vpcConfig.subnetwork.flatMap { s =>
          projectLabels.labels.collectFirst {
            case (labelName, labelValue) if labelName.equals(s) => labelValue
          }
        }

        VpcAndSubnetworkProjectLabelValues(vpcNetworkLabelValue, subnetworkLabelOption)
      }

    def fetchVpcLabelsFromProjectMetadata(
      vpcConfig: VirtualPrivateCloudLabels
    ): Future[Option[VpcAndSubnetworkProjectLabelValues]] =
      for {
        projectMetadataResponse <- projectMetadataRequest(vpcConfig).map(_.executeAsync().get())
        projectLabels <- projectMetadataResponseToLabels(projectMetadataResponse)
      } yield networkLabelsFromProjectLabels(vpcConfig, projectLabels)

    /*
    First, try to fetch the network information from labels, where that fetch may still return None.
    Then, if we did not discover a network via labels for whatever reason try to look for literal values.
     */
    def fetchVpcLabels(
      vpcConfig: VirtualPrivateCloudConfiguration
    ): Future[Option[VpcAndSubnetworkProjectLabelValues]] = {
      // Added explicit types to hopefully help future devs who stumble across this two-step code
      val fetchedFromLabels: Future[Option[VpcAndSubnetworkProjectLabelValues]] = vpcConfig.labelsOption match {
        case Some(labels: VirtualPrivateCloudLabels) => fetchVpcLabelsFromProjectMetadata(labels)
        case None => Future.successful(None)
      }
      fetchedFromLabels map {
        _ orElse {
          vpcConfig.literalsOption map { literals: VirtualPrivateCloudLiterals =>
            VpcAndSubnetworkProjectLabelValues(literals.network, literals.subnetwork)
          }
        }
      }
    }

    val vpcConfig: VirtualPrivateCloudConfiguration =
      gcpBatchConfiguration.batchAttributes.virtualPrivateCloudConfiguration
    fetchVpcLabels(vpcConfig)
  }

  override lazy val workflowPaths: Future[GcpBatchWorkflowPaths] = for {
    gcsCred <- gcsCredentials
    genomicsCred <- genomicsCredentials
    validatedPathBuilders <- pathBuilders
  } yield new GcpBatchWorkflowPaths(workflowDescriptor,
                                    gcsCred,
                                    genomicsCred,
                                    gcpBatchConfiguration,
                                    validatedPathBuilders,
                                    standardStreamNameToFileNameMetadataMapper
  )(ioEc)

  override lazy val initializationData: Future[GcpBackendInitializationData] = for {
    batchWorkflowPaths <- workflowPaths
    gcsCreds <- gcsCredentials
    vpcNetworkAndSubnetworkProjectLabels <- vpcNetworkAndSubnetworkProjectLabelsFuture()
  } yield models.GcpBackendInitializationData(
    workflowPaths = batchWorkflowPaths,
    runtimeAttributesBuilder = runtimeAttributesBuilder,
    gcpBatchConfiguration = gcpBatchConfiguration,
    gcsCredentials = gcsCreds,
    requestFactory = new GcpBatchRequestFactoryImpl()(gcpBatchConfiguration.batchAttributes.gcsTransferConfiguration),
    privateDockerEncryptionKeyName = privateDockerEncryptionKeyName,
    privateDockerEncryptedToken = privateDockerEncryptedToken,
    vpcNetworkAndSubnetworkProjectLabels = vpcNetworkAndSubnetworkProjectLabels
  )

  override def validateWorkflowOptions(): Try[Unit] = GcpLabel.fromWorkflowOptions(workflowOptions).map(_ => ())

  override def beforeAll(): Future[Option[BackendInitializationData]] =
    for {
      paths <- workflowPaths
      _ = publishWorkflowRoot(paths.workflowRoot.pathAsString)
      data <- initializationData
    } yield Option(data)

  def standardStreamNameToFileNameMetadataMapper(gcpBatchJobPaths: GcpBatchJobPaths, streamName: String): String =
    GcpBatchInitializationActor.defaultStandardStreamNameToFileNameMetadataMapper(gcpBatchJobPaths, streamName)

  override lazy val ioCommandBuilder: GcsBatchCommandBuilder.type = GcsBatchCommandBuilder

}

object GcpBatchInitializationActor {
  // For metadata publishing purposes default to using the name of a standard stream as the stream's filename.
  def defaultStandardStreamNameToFileNameMetadataMapper(gcpBatchJobPaths: GcpBatchJobPaths,
                                                        streamName: String
  ): String = streamName

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
