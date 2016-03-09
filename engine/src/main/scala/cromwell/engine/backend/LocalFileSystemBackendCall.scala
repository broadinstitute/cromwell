package cromwell.engine.backend

import java.nio.file.Path

/**
 * This represents any BackendCall that has access to the same file system as Cromwell itself
 */
trait LocalFileSystemBackendCall extends BackendCall {
  val returnCode: Path
}
