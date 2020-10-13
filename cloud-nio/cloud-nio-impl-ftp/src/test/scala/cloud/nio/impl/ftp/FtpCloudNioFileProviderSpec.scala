package cloud.nio.impl.ftp

import java.nio.channels.Channels

import cloud.nio.util.TryWithResource._
import common.assertion.CromwellTimeoutSpec
import org.mockftpserver.fake.filesystem.{DirectoryEntry, FileEntry}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class FtpCloudNioFileProviderSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with MockFtpFileSystem {

  behavior of "FtpCloudNioFileProviderSpec"

  lazy val fileProvider = mockProvider.fileProvider

  it should "check file existence" in {
    val root = "/root/existsPath/"
    val file = s"$root/file"
    fakeUnixFileSystem.add(new DirectoryEntry(root))
    fileProvider.existsPath("localhost", file) shouldBe false
    fakeUnixFileSystem.add(new FileEntry(file))
    fileProvider.existsPath("localhost", file) shouldBe true
  }

  it should "check directory existence" in {
    val root = "/root/existsPaths"
    val file = s"$root/file"
    val directory = s"$root/directory"
    fakeUnixFileSystem.add(new DirectoryEntry(root))
    fileProvider.existsPaths("localhost", directory) shouldBe false
    // Add a file, should still return false
    fakeUnixFileSystem.add(new FileEntry(file))
    fileProvider.existsPaths("localhost", file) shouldBe false
    // Now add the directory
    fakeUnixFileSystem.add(new DirectoryEntry(directory))
    fileProvider.existsPaths("localhost", directory) shouldBe true
  }

  it should "list objects" in {
    val root = "/root/listObjects"
    val file = s"$root/file"
    val directory = s"$root/directory"
    fakeUnixFileSystem.add(new DirectoryEntry(root))
    fakeUnixFileSystem.add(new FileEntry(file))
    fakeUnixFileSystem.add(new DirectoryEntry(directory))

    fileProvider.listObjects("localhost", root, None).paths should contain theSameElementsAs List(
      file.stripPrefix("/"), directory.stripPrefix("/")
    )
  }

  it should "copy files" in {
    val root = "/root/copy"
    val fileA = s"$root/fileA"
    val fileB = s"$root/fileB"
    fakeUnixFileSystem.add(new DirectoryEntry(root))
    val fileEntryA = new FileEntry(fileA)
    fileEntryA.setContents("salut!")
    fakeUnixFileSystem.add(fileEntryA)

    fileProvider.copy("localhost", fileA, "localhost", fileB)

    fakeUnixFileSystem.exists(fileB) shouldBe true
    Option(fakeUnixFileSystem.getEntry(fileB)) match {
      case Some(f: FileEntry) => f.getSize shouldBe 6
      case _ => fail("Copy failed")
    }
  }
  
  it should "throw when trying to copy files across hosts" in {
    an[UnsupportedOperationException] shouldBe thrownBy(fileProvider.copy("a", "b", "c", "d"))
  }

  it should "delete files" in {
    val root = "/root/deleteIfExists"
    val fileA = s"$root/fileA"
    val fileB = s"$root/fileB"
    fakeUnixFileSystem.add(new DirectoryEntry(root))
    fakeUnixFileSystem.add(new FileEntry(fileA))

    fileProvider.deleteIfExists("localhost", fileA) shouldBe true
    fileProvider.deleteIfExists("localhost", fileB) shouldBe false
  }

  it should "read a file" in {
    val root = "/root/read"
    val file = s"$root/file"
    fakeUnixFileSystem.add(new DirectoryEntry(root))
    val fileEntry = new FileEntry(file)
    fileEntry.setContents("salut!")
    fakeUnixFileSystem.add(fileEntry)

    val charBuffer = new Array[Char](6)
    val readableByteChannel = fileProvider.read("localhost", file, 0L)
    val reader = Channels.newReader(readableByteChannel, "ASCII")
    reader.read(charBuffer)
    charBuffer.mkString shouldBe "salut!"
    reader.close()
  }

  it should "write to a file" in {
    val root = "/root/write"
    val file = s"$root/file"
    fakeUnixFileSystem.add(new DirectoryEntry(root))

    val writableByteChannel = fileProvider.write("localhost", file)
    val writer = Channels.newOutputStream(writableByteChannel)
    writer.write("salut!".getBytes("ASCII"))
    writer.close()
    
    Option(fakeUnixFileSystem.getEntry(file)) match {
      case Some(f: FileEntry) =>
        var content: String = ""
        tryWithResource(() => f.createInputStream()) { is =>
          val byteBuffer = new Array[Byte](6)
          is.read(byteBuffer)
          content = byteBuffer.map(_.toChar).mkString
        }

        content shouldBe "salut!"
      case _ => fail("Failed to write")
    }
  }

  it should "get file attributes" in {
    val root = "/root/fileAttributes"
    val file = s"$root/file"
    fakeUnixFileSystem.add(new DirectoryEntry(root))
    val fileEntry = new FileEntry(file)
    fileEntry.setContents("salut!")
    fakeUnixFileSystem.add(fileEntry)

    val maybeAttributes = fileProvider.fileAttributes("localhost", file)
    maybeAttributes.isDefined shouldBe true
    val attributes = maybeAttributes.get
    attributes.size() shouldBe 6L
    attributes.isDirectory shouldBe false
    attributes.isRegularFile shouldBe true
    attributes.fileHash shouldBe None
    attributes.fileKey shouldBe "localhost/root/fileAttributes/file"
    attributes.lastModifiedTime() should not be null
  }

  it should "create a directory" in {
    val root = "/root/directory"
    val directory = s"$root/directory"
    fakeUnixFileSystem.add(new DirectoryEntry(root))
    fakeUnixFileSystem.exists(directory) shouldBe false
    fileProvider.createDirectory("localhost", directory)
    fakeUnixFileSystem.exists(directory) shouldBe true
  }
}
