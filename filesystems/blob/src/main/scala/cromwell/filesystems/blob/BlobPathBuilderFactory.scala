package cromwell.filesystems.blob

import akka.actor.ActorSystem
import cats.syntax.validated._
import com.azure.core.credential.TokenRequestContext
import com.azure.core.management.AzureEnvironment
import com.azure.core.management.profile.AzureProfile
import com.azure.identity.DefaultAzureCredentialBuilder
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import cromwell.core.path.PathBuilderFactory.PriorityBlob
import net.ceedubs.ficus.Ficus._
import net.ceedubs.ficus.readers.ArbitraryTypeReader._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.jdk.DurationConverters._
import scala.util.{Failure, Success, Try}


// WSM config is needed for accessing WSM-managed blob containers created in Terra workspaces.
// If the identity executing Cromwell has native access to the blob container, this can be ignored.
final case class WorkspaceManagerConfig(
    url: WorkspaceManagerURL,
    workspaceId: WorkspaceId,
    containerResourceId: ContainerResourceId,
    b2cToken: Option[String] // dev-only
)

final case class SubscriptionId(value: String) {override def toString: String = value}
final case class BlobContainerName(value: String) {override def toString: String = value}
final case class StorageAccountName(value: String) {override def toString: String = value}
final case class EndpointURL(value: String) {override def toString: String = value}
final case class WorkspaceId(value: String) {override def toString: String = value}
final case class ContainerResourceId(value: String) {override def toString: String = value}
final case class WorkspaceManagerURL(value: String) {override def toString: String = value}

final case class BlobPathBuilderFactory(globalConfig: Config, instanceConfig: Config, singletonConfig: Config) extends PathBuilderFactory {
  val subscription: Option[SubscriptionId] = instanceConfig.as[Option[String]]("subscription").map(SubscriptionId)
  val container: BlobContainerName = BlobContainerName(instanceConfig.as[String]("container"))
  val endpoint: EndpointURL = EndpointURL(instanceConfig.as[String]("endpoint"))
  val expiryBufferMinutes: Long = instanceConfig.as[Option[Long]]("expiry-buffer-minutes").getOrElse(10)
  val workspaceManagerConfig: Option[WorkspaceManagerConfig] = instanceConfig.as[Option[WorkspaceManagerConfig]]("workspace-manager")

  val blobTokenGenerator: BlobTokenGenerator = workspaceManagerConfig.map { wsmConfig =>
    val wsmClient: WorkspaceManagerApiClientProvider = new HttpWorkspaceManagerClientProvider(wsmConfig.url, wsmConfig.b2cToken.get)
    // parameterizing client instead of URL to make injecting mock client possible
    BlobTokenGenerator.createBlobTokenGenerator(container, endpoint, wsmConfig.workspaceId, wsmConfig.containerResourceId, wsmClient)
  }.getOrElse(
    BlobTokenGenerator.createBlobTokenGenerator(container, endpoint, subscription)
  )

  val fsm: BlobFileSystemManager = BlobFileSystemManager(container, endpoint, expiryBufferMinutes, blobTokenGenerator)

  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[BlobPathBuilder] = {
    Future {
      new BlobPathBuilder(container, endpoint)(fsm)
    }
  }

  override def priority: Int = PriorityBlob
}

/**
  * Strategy for obtaining an access token in an environment with available Azure identity.
  * If you need to disambiguate among multiple active user-assigned managed identities, pass
  * in the client id of the identity that should be used.
  */
case class AzureCredentials(identityClientId: Option[String]) {

  final val tokenAcquisitionTimeout = 30.seconds

  val azureProfile = new AzureProfile(AzureEnvironment.AZURE)
  val tokenScope = "https://management.azure.com/.default"

  def tokenRequestContext: TokenRequestContext = {
    val trc = new TokenRequestContext()
    trc.addScopes(tokenScope)
    trc
  }

  def defaultCredentialBuilder: DefaultAzureCredentialBuilder =
    new DefaultAzureCredentialBuilder()
      .authorityHost(azureProfile.getEnvironment.getActiveDirectoryEndpoint)

  def getAccessToken: ErrorOr[String] = {
    val credentials = identityClientId.foldLeft(defaultCredentialBuilder) {
      (builder, clientId) => builder.managedIdentityClientId(clientId)
    }.build()

    Try(
      credentials
        .getToken(tokenRequestContext)
        .block(tokenAcquisitionTimeout.toJava)
    ) match {
      case Success(null) => "null token value attempting to obtain access token".invalidNel
      case Success(token) => token.getToken.validNel
      case Failure(error) => s"Failed to refresh access token: ${error.getMessage}".invalidNel
    }
  }
}
