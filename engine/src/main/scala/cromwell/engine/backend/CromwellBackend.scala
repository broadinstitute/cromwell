package cromwell.engine.backend

import akka.actor.ActorSystem
import cromwell.core.WorkflowOptions

import scala.language.postfixOps
import scala.util.Try

/**
  * Provides a global singleton access to the instantiated backend. This isn't bulletproof and can lead to
  * exceptions being thrown but the way the system utilizes it it shouldn't be a problem (famous last words).
  *
  * The singleton access isn't a great thing long term however since we know that in the not too distant
  * future Cromwell will be sitting on multiple backends at once we'll need to redesign the innards to not
  * expect a single backend anyways. In terms of immediate-term grossness this was the least ugly of the possibilities
  * for now
  */
object CromwellBackend {
  private var _backends: Option[Map[String, Backend]] = None
  private var _defaultBackend: Option[Backend] = None

  def backendFromConfig(entry: BackendConfigurationEntry, actorSystem: ActorSystem): Backend = {
    val constructor = Class.forName(entry.className).getConstructor(classOf[BackendConfigurationEntry], classOf[ActorSystem])
    constructor.newInstance(entry, actorSystem).asInstanceOf[Backend]
  }

  def initBackends(backendEntries: List[BackendConfigurationEntry], defaultBackendEntry: BackendConfigurationEntry, actorSystem: ActorSystem) = {
    case class EntryAndBackend(entry: BackendConfigurationEntry, backend: Backend)
    val entriesAndBackends = backendEntries map { e => EntryAndBackend(e, backendFromConfig(e, actorSystem)) }
    _defaultBackend = entriesAndBackends find { _.entry == defaultBackendEntry } map { _.backend }
    _backends = Option(entriesAndBackends map { eb => eb.entry.name -> eb.backend } toMap)
  }

  def backend(backendName: String): Backend = _backends map { _(backendName) } getOrElse { throw new IllegalStateException("backend() called prior to initBackends") }

  def tryBackend(backendName: String) = Try(backend(backendName))

  /**
    * Register a new custom backend in the CromwellBackend set. This is mainly useful just for test cases which want
    * to use a custom/mock/stub backend in a test.
    *
    */
  private[cromwell] def registerCustomBackend(backendName: String, backend: Backend) = {
    _backends match {
      case Some(backends) => _backends = Some(backends + (backendName -> backend))
      case None => throw new IllegalStateException("registerCustomBackend(name, backend) called prior to initBackends")
    }
  }

  def defaultBackend = _defaultBackend getOrElse { throw new IllegalStateException("defaultBackend() called prior to initBackends") }

  // PBE:This is shorthand for reading a JSON options file and returning a backend.
  // Probably won't make it past the imminent move from workflow-scoped backend to call-scoped backend.
  def getBackendFromOptions(optionsString: String): Backend = {
    val triedBackend = for {
      workflowOptions <- WorkflowOptions.fromJsonString(optionsString)
      backendRequested <- workflowOptions.get("backend")
      backend <- tryBackend(backendRequested)
    } yield backend

    triedBackend getOrElse CromwellBackend.defaultBackend
  }
}
