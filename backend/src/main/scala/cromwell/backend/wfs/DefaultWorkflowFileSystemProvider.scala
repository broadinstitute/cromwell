package cromwell.backend.wfs

import java.nio.file.FileSystems

object DefaultWorkflowFileSystemProvider extends WorkflowFileSystemProvider {
  override def fileSystemOption(params: WorkflowFileSystemProviderParams) = {
    Option(FileSystems.getDefault)
  }
}
