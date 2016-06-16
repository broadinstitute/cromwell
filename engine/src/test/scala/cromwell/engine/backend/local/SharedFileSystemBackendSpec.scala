package cromwell.engine.backend.local

import java.nio.file.{FileSystems, Files, NoSuchFileException}

import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.CromwellSpec.IntegrationTest
import cromwell.filesystems.gcs.GoogleAuthMode.GoogleAuthOptions
import cromwell.filesystems.gcs._
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito

import scala.util.{Failure, Success}

class SharedFileSystemBackendSpec extends FlatSpec with Matchers with Mockito {

  behavior of "SharedFileSystemBackend"

  it should "localize the same path as already localized" ignore {
//    val orig = File.newTemp("file").write("hello world")
//    val workflowDescriptor = mock[OldStyleWorkflowDescriptor]
//
//    val result = SharedFileSystemBackend.localizePathAlreadyLocalized(orig.fullPath, orig.path, workflowDescriptor)
//    result should be(Success(()))
//
//    countLinks(orig) should be(1)
//
//    orig.delete(ignoreIOExceptions = true)
  }

  it should "localize a path already localized" ignore {
//    val orig = File.newTemp("file").write("hello world")
//    val dest = File.newTemp("file").delete(ignoreIOExceptions = false)
//    val workflowDescriptor = mock[OldStyleWorkflowDescriptor]
//
//    dest.toJava shouldNot exist
//    orig.copyTo(dest)
//    dest.toJava should exist
//    val result = SharedFileSystemBackend.localizePathAlreadyLocalized(orig.fullPath, dest.path, workflowDescriptor)
//    result should be(Success(()))
//
//    countLinks(orig) should be(1)
//    countLinks(dest) should be(1)
//
//    orig.delete(ignoreIOExceptions = true)
//    dest.delete(ignoreIOExceptions = true)
  }

  it should "not localize a path not localized" ignore {
//    val orig = File.newTemp("file").write("hello world")
//    val dest = File.newTemp("file").delete(ignoreIOExceptions = false)
//    val workflowDescriptor = mock[OldStyleWorkflowDescriptor]
//
//    dest.toJava shouldNot exist
//    val result = SharedFileSystemBackend.localizePathAlreadyLocalized(orig.fullPath, dest.path, workflowDescriptor)
//    result.isFailure should be(true)
//
//    val exception = intercept[NoSuchFileException](result.get)
//    exception.getMessage should endWith("does not exist")
//
//    orig.delete(ignoreIOExceptions = true)
//    dest.delete(ignoreIOExceptions = true)
  }

  it should "localize a path via copy" in {
//    val orig = File.newTemp("file").write("hello world")
//    val dest = File.newTemp("file").delete(ignoreIOExceptions = false)
//    val workflowDescriptor = mock[OldStyleWorkflowDescriptor]
//
//    dest.toJava shouldNot exist
//    val result = SharedFileSystemBackend.localizePathViaCopy(orig.fullPath, dest.path, workflowDescriptor)
//    result should be(Success(()))
//
//    dest.append(".")
//
//    orig.contentAsString should be("hello world")
//    countLinks(orig) should be(1)
//    countLinks(dest) should be(1)
//    isSymLink(dest) should be(false)
//    dest.contentAsString should be("hello world.")
//
//    orig.delete(ignoreIOExceptions = true)
//    dest.delete(ignoreIOExceptions = true)
  }

  it should "localize a path via hard link" ignore {
//    val orig = File.newTemp("file").write("hello world")
//    val dest = File.newTemp("file").delete(ignoreIOExceptions = false)
//    val workflowDescriptor = mock[OldStyleWorkflowDescriptor]
//
//    dest.toJava shouldNot exist
//    val result = SharedFileSystemBackend.localizePathViaHardLink(orig.fullPath, dest.path, workflowDescriptor)
//    result should be(Success(()))
//
//    dest.append(".")
//
//    orig.contentAsString should be("hello world.")
//    countLinks(orig) should be(2)
//    countLinks(dest) should be(2)
//    isSymLink(dest) should be(false)
//    dest.contentAsString should be("hello world.")
//
//    orig.delete(ignoreIOExceptions = true)
//    dest.delete(ignoreIOExceptions = true)
  }

  it should "localize a path via symbolic link" ignore {
//    val orig = File.newTemp("file").write("hello world")
//    val dest = File.newTemp("file").delete(ignoreIOExceptions = false)
//    val workflowDescriptor = mock[OldStyleWorkflowDescriptor]
//
//    dest.toJava shouldNot exist
//    val result = SharedFileSystemBackend.localizePathViaSymbolicLink(orig.fullPath, dest.path, workflowDescriptor)
//    result should be(Success(()))
//
//    dest.append(".")
//
//    orig.contentAsString should be("hello world.")
//    countLinks(orig) should be(1)
//    countLinks(dest) should be(1)
//    isSymLink(dest) should be(true)
//    dest.contentAsString should be("hello world.")
//
//    orig.delete(ignoreIOExceptions = true)
//    dest.delete(ignoreIOExceptions = true)
  }

  it should "localize from gcs" taggedAs IntegrationTest ignore {
    // via https://cloud.google.com/storage/docs/access-public-data
//    val origPath = "gs://uspto-pair/applications/05900002.zip"
//    val dest = File.newTemp("file").delete(ignoreIOExceptions = false)
//    val workflowDescriptor = mock[OldStyleWorkflowDescriptor]
//
//    val auth: GoogleAuthMode = new ApplicationDefaultMode("default")
//    val storage = auth.buildStorage(new GoogleAuthOptions {
//      override def get(key: String) = Failure(new UnsupportedOperationException("empty options"))
//    }, ConfigFactory.load)
//    val fileSystem = GcsFileSystemProvider(storage).getFileSystem
//    workflowDescriptor.fileSystems returns List(fileSystem, FileSystems.getDefault)
//
//    dest.toJava shouldNot exist
//    val result = SharedFileSystemBackend.localizeFromGcs(origPath, dest.path, workflowDescriptor)
//    result should be(Success(()))
//
//    countLinks(dest) should be(1)
//    isSymLink(dest) should be(false)
//    dest.toJava should have length 1811
//
//    dest.delete(ignoreIOExceptions = true)
  }

  it should "localize the same path as already localized in a directory" ignore {
//    val origDir = File.newTempDir("orig")
//    val orig = origDir.createChild("subdir/file").write("hello world")
//    val workflowDescriptor = mock[OldStyleWorkflowDescriptor]
//
//    val result = SharedFileSystemBackend.localizePathAlreadyLocalized(orig.fullPath, orig.path, workflowDescriptor)
//    result should be(Success(()))
//
//    countLinks(orig) should be(1)
//
//    origDir.delete(ignoreIOExceptions = true)
  }

  it should "localize a path already localized in a directory" ignore {
//    val origDir = File.newTempDir("orig")
//    val orig = origDir.createChild("subdir/file").write("hello world")
//    val destDir = File.newTempDir("dest")
//    val dest = destDir./("subdir/file")
//    destDir.delete(ignoreIOExceptions = false)
//    val workflowDescriptor = mock[OldStyleWorkflowDescriptor]
//
//    destDir.toJava shouldNot exist
//    dest.toJava shouldNot exist
//    dest.parent.createDirectories()
//    orig.copyTo(dest)
//    destDir.toJava should exist
//    dest.toJava should exist
//    val result = SharedFileSystemBackend.localizePathAlreadyLocalized(orig.fullPath, dest.path, workflowDescriptor)
//    result should be(Success(()))
//
//    countLinks(orig) should be(1)
//    countLinks(dest) should be(1)
//
//    origDir.delete(ignoreIOExceptions = true)
//    destDir.delete(ignoreIOExceptions = true)
  }

  it should "not localize a path not localized in a directory" ignore {
//    val origDir = File.newTempDir("orig")
//    val orig = origDir.createChild("subdir/file").write("hello world")
//    val destDir = File.newTempDir("dest")
//    val dest = destDir./("subdir/file")
//    destDir.delete(ignoreIOExceptions = false)
//    destDir.toJava shouldNot exist
//    dest.toJava shouldNot exist
//    val workflowDescriptor = mock[OldStyleWorkflowDescriptor]
//    val result = SharedFileSystemBackend.localizePathAlreadyLocalized(orig.fullPath, dest.path, workflowDescriptor)
//    result.isFailure should be(true)
//
//    val exception = intercept[NoSuchFileException](result.get)
//    exception.getMessage should endWith("does not exist")
//
//    origDir.delete(ignoreIOExceptions = true)
//    destDir.delete(ignoreIOExceptions = true)
  }

  it should "localize a path via copy in a directory" ignore {
//    val origDir = File.newTempDir("orig")
//    val orig = origDir.createChild("subdir/file").write("hello world")
//    val destDir = File.newTempDir("dest")
//    val dest = destDir./("subdir/file")
//    destDir.delete(ignoreIOExceptions = false)
//    val workflowDescriptor = mock[OldStyleWorkflowDescriptor]
//
//    destDir.toJava shouldNot exist
//    dest.toJava shouldNot exist
//    val result = SharedFileSystemBackend.localizePathViaCopy(orig.fullPath, dest.path, workflowDescriptor)
//    result should be(Success(()))
//
//    dest.append(".")
//
//    orig.contentAsString should be("hello world")
//    countLinks(orig) should be(1)
//    countLinks(dest) should be(1)
//    isSymLink(dest) should be(false)
//    dest.contentAsString should be("hello world.")
//
//    origDir.delete(ignoreIOExceptions = true)
//    destDir.delete(ignoreIOExceptions = true)
  }

  it should "localize a path via hard link in a directory" ignore {
//    val origDir = File.newTempDir("orig")
//    val orig = origDir.createChild("subdir/file").write("hello world")
//    val destDir = File.newTempDir("dest")
//    val dest = destDir./("subdir/file")
//    destDir.delete(ignoreIOExceptions = false)
//    val workflowDescriptor = mock[OldStyleWorkflowDescriptor]
//
//    destDir.toJava shouldNot exist
//    dest.toJava shouldNot exist
//    val result = SharedFileSystemBackend.localizePathViaHardLink(orig.fullPath, dest.path, workflowDescriptor)
//    result should be(Success(()))
//
//    dest.append(".")
//
//    orig.contentAsString should be("hello world.")
//    countLinks(orig) should be(2)
//    countLinks(dest) should be(2)
//    isSymLink(dest) should be(false)
//    dest.contentAsString should be("hello world.")
//
//    origDir.delete(ignoreIOExceptions = true)
//    destDir.delete(ignoreIOExceptions = true)
  }

  it should "localize a path via symbolic link in a directory" ignore {
//    val origDir = File.newTempDir("orig")
//    val orig = origDir.createChild("subdir/file").write("hello world")
//    val destDir = File.newTempDir("dest")
//    val dest = destDir./("subdir/file")
//    destDir.delete(ignoreIOExceptions = false)
//    val workflowDescriptor = mock[OldStyleWorkflowDescriptor]
//
//    destDir.toJava shouldNot exist
//    dest.toJava shouldNot exist
//    val result = SharedFileSystemBackend.localizePathViaSymbolicLink(orig.fullPath, dest.path, workflowDescriptor)
//    result should be(Success(()))
//
//    dest.append(".")
//
//    orig.contentAsString should be("hello world.")
//    countLinks(orig) should be(1)
//    countLinks(dest) should be(1)
//    isSymLink(dest) should be(true)
//    dest.contentAsString should be("hello world.")
//
//    origDir.delete(ignoreIOExceptions = true)
//    destDir.delete(ignoreIOExceptions = true)
  }

  it should "localize from gcs to a directory" taggedAs IntegrationTest ignore {
    // via https://cloud.google.com/storage/docs/access-public-data
//    val origPath = "gs://uspto-pair/applications/05900002.zip"
//    val destDir = File.newTempDir("dest")
//    val dest = destDir./("subdir/file")
//    destDir.delete(ignoreIOExceptions = false)
//    val auth: GoogleAuthMode = new ApplicationDefaultMode("default")
//    val storage = auth.buildStorage(new GoogleAuthOptions {
//      override def get(key: String) = Failure(new UnsupportedOperationException("empty options"))
//    }, ConfigFactory.load())
//    val fileSystem = GcsFileSystemProvider(storage).getFileSystem
//    val workflowDescriptor = mock[OldStyleWorkflowDescriptor]
//    workflowDescriptor.fileSystems returns List(fileSystem, FileSystems.getDefault)
//
//    destDir.toJava shouldNot exist
//    dest.toJava shouldNot exist
//    val result = SharedFileSystemBackend.localizeFromGcs(origPath, dest.path, workflowDescriptor)
//    result should be(Success(()))
//
//    countLinks(dest) should be(1)
//    isSymLink(dest) should be(false)
//    dest.toJava should have length 1811
//
//    destDir.delete(ignoreIOExceptions = true)
  }

  private[this] def countLinks(file: File): Int = Files.getAttribute(file.path, "unix:nlink").asInstanceOf[Int]

  private[this] def isSymLink(file: File): Boolean = Files.isSymbolicLink(file.path)
}
