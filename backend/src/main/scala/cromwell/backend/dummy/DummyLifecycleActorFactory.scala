package cromwell.backend.dummy

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.standard.callcaching.{StandardCacheHitCopyingActor, StandardFileHashingActor}
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardInitializationActor, StandardLifecycleActorFactory}

class DummyLifecycleActorFactory(override val name: String, override val configurationDescriptor: BackendConfigurationDescriptor) extends StandardLifecycleActorFactory {

  /**
    * @return the key to use for storing and looking up the job id.
    */
  override def jobIdKey: String = "__dummy_operation_id"

  /**
    * @return the asynchronous executor class.
    */
  override def asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] = classOf[DummyAsyncExecutionActor]

  // Don't cache-hit copy
  override lazy val cacheHitCopyingActorClassOption: Option[Class[_ <: StandardCacheHitCopyingActor]] = None

  // Don't hash files
  override lazy val fileHashingActorClassOption: Option[Class[_ <: StandardFileHashingActor]] = None

  override def backendSingletonActorProps(serviceRegistryActor: ActorRef): Option[Props] = Option(Props(new DummySingletonActor()))

  override lazy val initializationActorClass: Class[_ <: StandardInitializationActor] = classOf[DummyInitializationActor]

}
