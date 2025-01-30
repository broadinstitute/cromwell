package cromwell.backend.google.batch

import akka.actor.{ActorRef, Props}
import com.google.api.client.util.ExponentialBackOff
import com.google.api.gax.rpc.FixedHeaderProvider
import com.google.cloud.batch.v1.BatchServiceSettings
import com.google.common.collect.ImmutableMap
import com.typesafe.scalalogging.StrictLogging
import cromwell.backend._
import cromwell.backend.google.batch.GcpBatchBackendLifecycleActorFactory.preemptionCountKey
import cromwell.backend.google.batch.actors._
import cromwell.backend.google.batch.api.request.{BatchRequestExecutor, RequestHandler}
import cromwell.backend.google.batch.authentication.GcpBatchDockerCredentials
import cromwell.backend.google.batch.callcaching.{BatchBackendCacheHitCopyingActor, BatchBackendFileHashingActor}
import cromwell.backend.google.batch.models.{
  GcpBackendInitializationData,
  GcpBatchConfiguration,
  GcpBatchConfigurationAttributes
}
import cromwell.backend.standard._
import cromwell.backend.standard.callcaching.{StandardCacheHitCopyingActor, StandardFileHashingActor}
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.{CallOutputs, DockerCredentials}
import wom.graph.CommandCallNode

import scala.util.{Failure, Success, Try}

class GcpBatchBackendLifecycleActorFactory(override val name: String,
                                           override val configurationDescriptor: BackendConfigurationDescriptor
) extends StandardLifecycleActorFactory
    with GcpPlatform {

  override val requestedKeyValueStoreKeys: Seq[String] = Seq(preemptionCountKey)
  import GcpBatchBackendLifecycleActorFactory._

  override def jobIdKey: String = "__gcp_batch"
  protected val googleConfig: GoogleConfiguration = GoogleConfiguration(configurationDescriptor.globalConfig)

  override lazy val initializationActorClass: Class[_ <: StandardInitializationActor] =
    classOf[GcpBatchInitializationActor]

  override def asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[GcpBatchAsyncBackendJobExecutionActor]

  override lazy val finalizationActorClassOption: Option[Class[_ <: StandardFinalizationActor]] =
    Option(classOf[GcpBatchFinalizationActor])

  protected val batchAttributes: GcpBatchConfigurationAttributes = {
    def defaultBuildAttributes() =
      GcpBatchConfigurationAttributes(googleConfig, configurationDescriptor.backendConfig, "batchConfig")
    robustBuildAttributes(defaultBuildAttributes)
  }

  val batchConfiguration = new GcpBatchConfiguration(configurationDescriptor, googleConfig, batchAttributes)

  override def workflowInitializationActorParams(workflowDescriptor: BackendWorkflowDescriptor,
                                                 ioActor: ActorRef,
                                                 calls: Set[CommandCallNode],
                                                 serviceRegistryActor: ActorRef,
                                                 restart: Boolean
  ): StandardInitializationActorParams =
    GcpBatchInitializationActorParams(workflowDescriptor,
                                      ioActor,
                                      calls,
                                      batchConfiguration,
                                      serviceRegistryActor,
                                      restart
    )

  override def workflowFinalizationActorParams(workflowDescriptor: BackendWorkflowDescriptor,
                                               ioActor: ActorRef,
                                               calls: Set[CommandCallNode],
                                               jobExecutionMap: JobExecutionMap,
                                               workflowOutputs: CallOutputs,
                                               initializationDataOption: Option[BackendInitializationData]
  ): StandardFinalizationActorParams =
    GcpBatchFinalizationActorParams(workflowDescriptor,
                                    ioActor,
                                    batchConfiguration,
                                    calls,
                                    jobExecutionMap,
                                    workflowOutputs,
                                    initializationDataOption
    )

  override lazy val cacheHitCopyingActorClassOption: Option[Class[_ <: StandardCacheHitCopyingActor]] =
    Option(classOf[BatchBackendCacheHitCopyingActor])

  override lazy val fileHashingActorClassOption: Option[Class[_ <: StandardFileHashingActor]] = Option(
    classOf[BatchBackendFileHashingActor]
  )

  override def backendSingletonActorProps(serviceRegistryActor: ActorRef): Option[Props] = {
    implicit val requestHandler: RequestHandler = new RequestHandler

    val batchSettings = BatchServiceSettings.newBuilder
      .setHeaderProvider(FixedHeaderProvider.create(ImmutableMap.of("user-agent", "cromwell")))
      .build

    Option(
      GcpBatchBackendSingletonActor.props(
        qps = batchConfiguration.batchAttributes.qps,
        requestWorkers = batchConfiguration.batchAttributes.requestWorkers,
        serviceRegistryActor = serviceRegistryActor,
        batchRequestExecutor = new BatchRequestExecutor.CloudImpl(batchSettings)
      )(requestHandler)
    )
  }

  override def dockerHashCredentials(
    workflowDescriptor: BackendWorkflowDescriptor,
    initializationData: Option[BackendInitializationData]
  ): List[Any] =
    Try(BackendInitializationData.as[GcpBackendInitializationData](initializationData)) match {
      case Success(data) =>
        val tokenFromWorkflowOptions = workflowDescriptor.workflowOptions
          .get(GoogleAuthMode.DockerCredentialsTokenKey)
          .toOption
        val effectiveToken = tokenFromWorkflowOptions.orElse {
          data.gcpBatchConfiguration.dockerCredentials.map(_.token)
        }

        val dockerCredentials: Option[GcpBatchDockerCredentials] = effectiveToken.map { token =>
          // These credentials are being returned for hashing and all that matters in this context is the token
          // so just `None` the auth and key.
          val baseDockerCredentials = new DockerCredentials(token = token, authName = None, keyName = None)
          GcpBatchDockerCredentials.apply(baseDockerCredentials, googleConfig)
        }
        val googleCredentials = Option(data.gcsCredentials)
        List(dockerCredentials, googleCredentials).flatten

      case _ => List.empty[Any]
    }
}

object GcpBatchBackendLifecycleActorFactory extends StrictLogging {
  val preemptionCountKey = "PreemptionCount"

  private[batch] def robustBuildAttributes(buildAttributes: () => GcpBatchConfigurationAttributes,
                                           maxAttempts: Int = 3,
                                           initialIntervalMillis: Int = 5000,
                                           maxIntervalMillis: Int = 10000,
                                           multiplier: Double = 1.5,
                                           randomizationFactor: Double = 0.5
  ): GcpBatchConfigurationAttributes = {
    val backoff = new ExponentialBackOff.Builder()
      .setInitialIntervalMillis(initialIntervalMillis)
      .setMaxIntervalMillis(maxIntervalMillis)
      .setMultiplier(multiplier)
      .setRandomizationFactor(randomizationFactor)
      .build()

    // `attempt` is 1-based
    def build(attempt: Int): Try[GcpBatchConfigurationAttributes] =
      Try {
        buildAttributes()
      } recoverWith {
        // Try again if this was an Exception (as opposed to an Error) and we have not hit maxAttempts
        case ex: Exception if attempt < maxAttempts =>
          logger
            .warn(s"Failed to build GcpBatchConfigurationAttributes on attempt $attempt of $maxAttempts, retrying.", ex)
          Thread.sleep(backoff.nextBackOffMillis())
          build(attempt + 1)
        case e =>
          Failure(
            new RuntimeException(s"Failed to build GcpBatchConfigurationAttributes on attempt $attempt of $maxAttempts",
                                 e
            )
          )
      }
    // This intentionally throws if the final result of `build` is a `Failure`.
    build(attempt = 1).get
  }
}
