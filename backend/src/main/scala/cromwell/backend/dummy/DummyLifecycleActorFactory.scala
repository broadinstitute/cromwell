package cromwell.backend.dummy

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardLifecycleActorFactory}

class DummyLifecycleActorFactory(override val name: String, override val configurationDescriptor: BackendConfigurationDescriptor) extends StandardLifecycleActorFactory {

  /**
    * Returns the key to use for storing and looking up the job id.
    *
    * @return the key to use for storing and looking up the job id.
    */
  override def jobIdKey: String = "__dummy_operation_id"

  /**
    * Returns the asynchronous executor class.
    *
    * @return the asynchronous executor class.
    */
  override def asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] = classOf[DummyAsyncExecutionActor]

}
