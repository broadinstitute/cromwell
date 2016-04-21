package cromwell.engine.backend

import akka.actor.ActorSystem
import cromwell.backend.BackendLifecycleActorFactory
import cromwell.core.WorkflowOptions

import scala.language.postfixOps
import scala.util.{Success, Try, Failure}

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

  /** PBE: Once we can switch to only using the _shadow* vars instead of _backends and _defaultBackend,
    *      then we don't need to do initBackends() anymore, and the _shadow* vars can become vals
    */
  private var _backends: Option[Map[String, Backend]] = None
  private var _defaultBackend: Option[Backend] = None

  private var _shadowBackendFactories: Option[Map[String, BackendLifecycleActorFactory]] = None
  private var _shadowDefaultBackendName: String = ""

  private def backendFromConfig(entry: BackendConfigurationEntry, actorSystem: ActorSystem): Backend = {
    val constructor = Class.forName(entry.className).getConstructor(classOf[BackendConfigurationEntry], classOf[ActorSystem])
    constructor.newInstance(entry, actorSystem).asInstanceOf[Backend]
  }

  def initBackends(backendEntries: List[BackendConfigurationEntry],
                   defaultBackendEntry: BackendConfigurationEntry,
                   actorSystem: ActorSystem,
                   shadowExecutionEnabled: Boolean) = {
    case class EntryAndBackend(entry: BackendConfigurationEntry, backend: Backend)
    val entriesAndBackends = backendEntries map { e => EntryAndBackend(e, backendFromConfig(e, actorSystem)) }
    _defaultBackend = entriesAndBackends find { _.entry == defaultBackendEntry } map { _.backend }
    _backends = Option(entriesAndBackends map { eb => eb.entry.name -> eb.backend } toMap)

    if (shadowExecutionEnabled) {
      val backendLifecycleActorFactories = backendEntries.map(e => e.name -> e.asBackendLifecycleActorFactory).toMap
      _shadowBackendFactories = Option(backendLifecycleActorFactories)
      _shadowDefaultBackendName = defaultBackendEntry.name
    }
  }

  def shadowBackendLifecycleFactory(backendName: String): Try[BackendLifecycleActorFactory] = {
    _shadowBackendFactories.map(_.get(backendName)) match {
      case None => Failure(new IllegalStateException("initBackends() was not called"))
      case Some(None) => Failure(new Exception(s"Backend $backendName was not found"))
      case Some(Some(factory)) => Success(factory)
    }
  }

  def shadowDefaultBackend: String = _shadowDefaultBackendName

  def isValidBackendName(name: String): Boolean = _shadowBackendFactories.exists(_.contains(name))

  private def backend(backendName: String): Backend = _backends map { _(backendName) } getOrElse { throw new IllegalStateException("backend() called prior to initBackends") }

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

  private def defaultBackend = _defaultBackend getOrElse { throw new IllegalStateException("defaultBackend() called prior to initBackends") }

  // PBE:This is shorthand for reading a JSON options file and returning a backend.
  // Probably won't make it past the imminent move from workflow-scoped backend to call-scoped backend.
  def getBackendFromOptions(optionsString: String): Backend = {
    val triedBackend = for {
      workflowOptions <- WorkflowOptions.fromJsonString(optionsString)
      backendRequested <- workflowOptions.get("backend")
      backend <- Try(backend(backendRequested))
    } yield backend

    triedBackend getOrElse CromwellBackend.defaultBackend
  }
}
