package cromwell.engine.backend

import akka.actor.ActorSystem
import cromwell.backend.BackendLifecycleActorFactory

import scala.language.postfixOps
import scala.util.{Success, Try, Failure}

/**
  * Provides a global singleton access to the instantiated backend factories.
  */
case class CromwellBackends(backendEntries: List[BackendConfigurationEntry],
                            defaultBackendEntry: BackendConfigurationEntry,
                            actorSystem: ActorSystem) {

  val backendLifecycleActorFactories = backendEntries.map(e => e.name -> e.asBackendLifecycleActorFactory).toMap
  private val _shadowBackendFactories = backendLifecycleActorFactories
  private val _shadowDefaultBackendName = defaultBackendEntry.name

  def shadowBackendLifecycleFactory(backendName: String): Try[BackendLifecycleActorFactory] = {
    _shadowBackendFactories.get(backendName) match {
      case None => Failure(new Exception(s"Backend $backendName was not found"))
      case Some(factory) => Success(factory)
    }
  }

  def shadowDefaultBackend: String = _shadowDefaultBackendName

  def isValidBackendName(name: String): Boolean = _shadowBackendFactories.contains(name)
}

object CromwellBackends {

  var instance: Option[CromwellBackends] = None


  def isValidBackendName(name: String): Boolean = evaluateOrThrow(_.isValidBackendName(name))
  def shadowDefaultBackend = evaluateOrThrow(_.shadowDefaultBackend)
  def shadowBackendLifecycleFactory(backendName: String) = evaluateOrThrow(_.shadowBackendLifecycleFactory(backendName))

  private def evaluateOrThrow[A](func: CromwellBackends => A): A = {
    instance match {
      case Some(cromwellBackend) => func(cromwellBackend)
      case None => throw new Exception("Cannot use CromwellBackend until initBackends is called")
    }
  }

  def initBackends(backendEntries: List[BackendConfigurationEntry],
                   defaultBackendEntry: BackendConfigurationEntry,
                   actorSystem: ActorSystem): Unit = {
    instance = Option(CromwellBackends(
      backendEntries: List[BackendConfigurationEntry],
      defaultBackendEntry: BackendConfigurationEntry,
      actorSystem: ActorSystem))
  }
}
