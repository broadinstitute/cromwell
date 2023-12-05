package cromwell.engine.backend

import cats.syntax.option._
import common.util.TryUtil
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.BackendLifecycleActorFactory

/**
  * Provides a global singleton access to the instantiated backend factories.
  */
case class CromwellBackends(backendEntries: List[BackendConfigurationEntry]) {

  // Raise the exception here if some backend factories failed to instantiate
  val backendLifecycleActorFactories =
    TryUtil.sequenceMap(backendEntries.map(e => e.name -> e.asBackendLifecycleActorFactory).toMap).get

  def backendLifecycleActorFactoryByName(backendName: String): ErrorOr[BackendLifecycleActorFactory] =
    backendLifecycleActorFactories.get(backendName).toValidNel(s"Backend $backendName was not found")

  def isValidBackendName(name: String): Boolean = backendLifecycleActorFactories.contains(name)
}

object CromwellBackends {

  var instance: Option[CromwellBackends] = None

  def isValidBackendName(name: String): Boolean = evaluateIfInitialized(_.isValidBackendName(name))

  def backendLifecycleFactoryActorByName(backendName: String): ErrorOr[BackendLifecycleActorFactory] =
    evaluateIfInitialized(_.backendLifecycleActorFactoryByName(backendName))

  private def evaluateIfInitialized[A](func: CromwellBackends => A): A =
    instance match {
      case Some(cromwellBackend) => func(cromwellBackend)
      case None => throw new Exception("Cannot use CromwellBackend until initBackends is called")
    }

  def initBackends(backendEntries: List[BackendConfigurationEntry]): Unit =
    instance = Option(CromwellBackends(backendEntries))
}
