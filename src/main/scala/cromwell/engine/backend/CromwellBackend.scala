package cromwell.engine.backend

import akka.actor.ActorSystem
import org.slf4j.LoggerFactory
import Backend.BackendyString

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
  private val log = LoggerFactory.getLogger(getClass.getName)
  private var _backend: Option[Backend] = None

  def initBackend(backendType: String, actorSystem: ActorSystem): Backend = {
    _backend match {
      case None =>
        val backend = Backend.from(backendType, actorSystem)
        _backend = Option(backend)
        backend
      case Some(x) if x.backendType == backendType.toBackendType =>
        log.error(errorString(x, backendType))
        x
      case Some(x) => throw new IllegalStateException(errorString(x, backendType))
    }
  }

  def backend() = _backend match {
    case Some(x) => x
    case None => throw new IllegalStateException("Backend called prior to initBackend")
  }

  private def errorString(newBackend: Backend, oldBackend: String): String = {
    s"Backend already initialized to ${newBackend.backendType} attempting to change it to $oldBackend"
  }
}
