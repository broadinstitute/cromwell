package cromwell.filesystems.gcs

import cromwell.core.WorkflowOptions
import cromwell.filesystems.gcs.auth.GoogleAuthMode.NoAuthMode
import org.specs2.mock.Mockito

object MockGcsPathBuilder extends Mockito {
  val mockPathBuilder = GcsPathBuilderFactory(NoAuthMode).withOptions(mock[WorkflowOptions])
}
