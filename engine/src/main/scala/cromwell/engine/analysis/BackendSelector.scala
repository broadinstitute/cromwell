package cromwell.engine.analysis

import cromwell.engine.backend.Backend
import wdl4s.Call

import scala.util.Try

// This can be enhanced later...
object BackendSelector {
  def selectBackend(defaultBackend: Option[Backend], call: Call): Try[Backend] =
  // TODO: Probably want to enhance this somehow:
    Try(defaultBackend.get)
}
