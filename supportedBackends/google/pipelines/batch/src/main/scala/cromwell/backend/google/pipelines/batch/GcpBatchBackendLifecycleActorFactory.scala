package cromwell.backend.google.pipelines.batch

import akka.actor.{ActorRef, Props}
import com.google.api.client.util.ExponentialBackOff
import com.typesafe.scalalogging.StrictLogging
import cromwell.backend.BackendWorkflowDescriptor
import cromwell.backend.{BackendInitializationData, JobExecutionMap}
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.CallOutputs
import wom.graph.CommandCallNode
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.google.pipelines.batch.GcpBatchBackendLifecycleActorFactory.robustBuildAttributes
import cromwell.backend.standard._
import scala.util.{Failure, Try}


class GcpBatchBackendLifecycleActorFactory(name: String, override val configurationDescriptor: BackendConfigurationDescriptor)
  extends StandardLifecycleActorFactory {

  override def name: String = "batch"

  override def jobIdKey: String = "__gcp_batch"
  protected val googleConfig: GoogleConfiguration = GoogleConfiguration(configurationDescriptor.globalConfig)

  override lazy val initializationActorClass: Class[_ <: StandardInitializationActor] = classOf[GcpBatchInitializationActor]

  override def asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[GcpBatchAsyncBackendJobExecutionActor]

  override lazy val finalizationActorClassOption: Option[Class[_ <: StandardFinalizationActor]] =
    Option(classOf[GcpBatchFinalizationActor])



  protected val batchAttributes: GcpBatchConfigurationAttributes = {
    def defaultBuildAttributes()=
      GcpBatchConfigurationAttributes(googleConfig, configurationDescriptor.backendConfig, "batchConfig")
    robustBuildAttributes(defaultBuildAttributes)
  }

  val batchConfiguration = new GcpBatchConfiguration(configurationDescriptor, googleConfig, batchAttributes)

  override def workflowInitializationActorParams(

                                                  workflowDescriptor: BackendWorkflowDescriptor,
                                                  ioActor: ActorRef,
                                                  calls: Set[CommandCallNode],
                                                  serviceRegistryActor: ActorRef,
                                                  restart: Boolean): StandardInitializationActorParams = {
    GcpBatchInitializationActorParams(workflowDescriptor, ioActor , calls, batchConfiguration, serviceRegistryActor, restart)
  }

  override def workflowFinalizationActorParams(
                                           workflowDescriptor: BackendWorkflowDescriptor,
                                           ioActor: ActorRef,
                                           //batchConfiguration: GcpBatchConfiguration,
                                           calls: Set[CommandCallNode],
                                           jobExecutionMap: JobExecutionMap,
                                           workflowOutputs: CallOutputs,
                                           initializationDataOption: Option[BackendInitializationData]): StandardFinalizationActorParams = {
    GcpBatchFinalizationActorParams(workflowDescriptor, ioActor, batchConfiguration, calls, jobExecutionMap, workflowOutputs, initializationDataOption)
  }

  override def backendSingletonActorProps(serviceRegistryActor: ActorRef): Option[Props] = {
    Option(GcpBatchBackendSingletonActor.props("gcp-batch"))
  }

  //override def backendSingletonActorProps(serviceRegistryActor: ActorRef): Option[Props] = super
  //  .backendSingletonActorProps(serviceRegistryActor)

}

object GcpBatchBackendLifecycleActorFactory extends StrictLogging {

  private def robustBuildAttributes(buildAttributes: () => GcpBatchConfigurationAttributes,
                                            maxAttempts: Int = 3,
                                            initialIntervalMillis: Int = 5000,
                                            maxIntervalMillis: Int = 10000,
                                            multiplier: Double = 1.5,
                                            randomizationFactor: Double = 0.5): GcpBatchConfigurationAttributes = {
    val backoff = new ExponentialBackOff.Builder()
      .setInitialIntervalMillis(initialIntervalMillis)
      .setMaxIntervalMillis(maxIntervalMillis)
      .setMultiplier(multiplier)
      .setRandomizationFactor(randomizationFactor)
      .build()

    // `attempt` is 1-based
    def build(attempt: Int): Try[GcpBatchConfigurationAttributes] = {
      Try {
        buildAttributes()
      } recoverWith {
        // Try again if this was an Exception (as opposed to an Error) and we have not hit maxAttempts
        case ex: Exception if attempt < maxAttempts =>
          logger
            .warn(s"Failed to build GcpBatchConfigurationAttributes on attempt $attempt of $maxAttempts, retrying.", ex)
          Thread.sleep(backoff.nextBackOffMillis())
          build(attempt + 1)
        case e => Failure(new RuntimeException(s"Failed to build GcpBatchConfigurationAttributes on attempt $attempt of $maxAttempts", e))
      }
    }
    // This intentionally throws if the final result of `build` is a `Failure`.
    build(attempt = 1).get
  }


}
