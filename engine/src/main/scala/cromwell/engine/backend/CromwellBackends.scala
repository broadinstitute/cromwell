package cromwell.engine.backend

import cromwell.backend.BackendLifecycleActorFactory
import lenthall.util.TryUtil

import scala.util.{Failure, Success, Try}

/**
  * Provides a global singleton access to the instantiated backend factories.
  */
case class CromwellBackends(backendEntries: List[BackendConfigurationEntry]) {

  // Raise the exception here if some backend factories failed to instantiate
  val backendLifecycleActorFactories = TryUtil.sequenceMap(backendEntries.map(e => e.name -> e.asBackendLifecycleActorFactory).toMap).get

  def backendLifecycleActorFactoryByName(backendName: String): Try[BackendLifecycleActorFactory] = {
    backendLifecycleActorFactories.get(backendName) match {
      case None => Failure(new Exception(s"Backend $backendName was not found"))
      case Some(factory) => Success(factory)
    }
  }

  def isValidBackendName(name: String): Boolean = backendLifecycleActorFactories.contains(name)
}

object CromwellBackends {

  var instance: Option[CromwellBackends] = None


  def isValidBackendName(name: String): Boolean = evaluateIfInitialized(_.isValidBackendName(name))
  def backendLifecycleFactoryActorByName(backendName: String) = evaluateIfInitialized(_.backendLifecycleActorFactoryByName(backendName))

  private def evaluateIfInitialized[A](func: CromwellBackends => A): A = {
    instance match {
      case Some(cromwellBackend) => func(cromwellBackend)
      case None => throw new Exception("Cannot use CromwellBackend until initBackends is called")
    }
  }

  def initBackends(backendEntries: List[BackendConfigurationEntry]): Unit = {
    instance = Option(CromwellBackends(backendEntries))
  }
}
