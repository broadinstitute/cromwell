package cromwell.backend.impl.jes

import java.nio.file.Path

import cromwell.filesystems.gcs.GcsFileSystem
import org.mockito
import org.specs2.mock.Mockito

trait MockGcsFileSystemBuilder { self: Mockito =>

  protected def buildMockedGcsFileSystem = {
    val mockedPath = mock[Path]
    mockedPath.resolve(anyString) returns mockedPath
    val mockedGcsFileSystem = mock[GcsFileSystem]
    mockedGcsFileSystem.getPath(anyString, mockito.Matchers.anyVararg[String]()) returns mockedPath
    mockedGcsFileSystem
  }

}
