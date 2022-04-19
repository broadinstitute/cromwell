package cromwell.filesystems.oss.nio

import java.net.URI
import java.nio.charset.Charset
import java.nio.file.{DirectoryStream, NoSuchFileException, Path, StandardOpenOption}

import cromwell.core.TestKitSuite
import org.scalatest.BeforeAndAfter

import scala.jdk.CollectionConverters._
import scala.collection.mutable.ArrayBuffer
import scala.util.control.Breaks

class OssStorageFileSystemProviderSpec extends TestKitSuite with OssNioUtilSpec with BeforeAndAfter {
  behavior of "OssStorageFileSystemProviderSpec"

  it should "has right schema" in {
    mockProvider.getScheme shouldEqual OssStorageFileSystem.URI_SCHEMA
  }

  it should "work when creating new file system" in {
    val fs = mockProvider.newFileSystem(URI.create(s"oss://$bucket"), mockOssConf.toMap.asJava)
    fs.bucket shouldEqual bucket

    an [IllegalArgumentException] should be thrownBy ossProvider.newFileSystem(URI.create(s"oss://"), mockOssConf.toMap.asJava)
    an [IllegalArgumentException] should be thrownBy ossProvider.newFileSystem(URI.create(s"oss://$bucket:8812"), mockOssConf.toMap.asJava)

    val fs1 = mockProvider.getFileSystem(URI.create(s"oss://$bucket"))
    fs1.bucket shouldEqual bucket
  }

  it should "work when getting a new oss path" in {
    val path = mockProvider.getPath(URI.create(s"oss://$bucket$fileName"))
    path.bucket shouldEqual bucket
    path.key shouldEqual fileName.stripPrefix(OssStorageFileSystem.SEPARATOR)
  }

  it should "work when creating an output stream" taggedAs NeedAK in {
    val path = ossProvider.getPath(URI.create(s"oss://$bucket$fileName"))

    val outS = ossProvider.newOutputStream(path)
    outS.write(fileContent.getBytes)

    contentAsString(path) shouldEqual fileContent
    outS.asInstanceOf[OssAppendOutputStream].position shouldEqual fileContent.length
  }

  it should "work when creating an byte channel" taggedAs NeedAK in {
    val path = ossProvider.getPath(URI.create(s"oss://$bucket$fileName"))

    val outS = ossProvider.newOutputStream(path)
    outS.write(fileContent.getBytes)

    val inC = ossProvider.newByteChannel(path, Set(StandardOpenOption.READ).asJava)

    import java.nio.ByteBuffer
    val buf = ByteBuffer.allocate(1)

    val loop = new Breaks
    val builder = new StringBuilder

    var bytesRead = inC.read(buf)
    loop.breakable {
      while (bytesRead != -1) {
        buf.flip()
        val charset = Charset.forName("UTF-8")

        builder.append(charset.decode(buf).toString)
        buf.clear
        bytesRead = inC.read(buf)
      }
    }

    builder.toString shouldEqual fileContent
  }

  it should "delete file if it exists" taggedAs NeedAK in {
    val path = ossProvider.getPath(URI.create(s"oss://$bucket$fileName"))

    val outS = ossProvider.newOutputStream(path)
    outS.write(fileContent.getBytes)
    outS.close()

    ossProvider.deleteIfExists(path) shouldEqual true
    ossProvider.deleteIfExists(path) shouldEqual false
    an [NoSuchFileException] should be thrownBy ossProvider.delete(path)
  }

  it should "work when copying an object" taggedAs NeedAK in {
    val src = ossProvider.getPath(URI.create(s"oss://$bucket$fileName"))
    val target = ossProvider.getPath(URI.create(s"oss://$bucket${fileName}1"))
    ossProvider.deleteIfExists(src)

    writeObject(src)

    ossProvider.copy(src, target)

    ossProvider.deleteIfExists(target) shouldEqual true
  }

  it should "work when moving an object" taggedAs NeedAK in {
    val src = ossProvider.getPath(URI.create(s"oss://$bucket$fileName"))
    val target = ossProvider.getPath(URI.create(s"oss://$bucket${fileName}1"))
    ossProvider.deleteIfExists(src)

    writeObject(src)

    ossProvider.move(src, target)

    ossProvider.deleteIfExists(target) shouldEqual true
    ossProvider.deleteIfExists(src) shouldEqual false

  }

  it should "work for some basic operations" taggedAs NeedAK in {
    val path = ossProvider.getPath(URI.create(s"oss://$bucket$fileName"))
    val path1 = ossProvider.getPath(URI.create(s"oss://$bucket$fileName"))

    ossProvider.isHidden(path) shouldEqual false
    ossProvider.isSameFile(path, path1)

    an [UnsupportedOperationException] should be thrownBy ossProvider.getFileStore(path)

    an [NoSuchFileException] should be thrownBy ossProvider.checkAccess(path)

    val dir = ossProvider.getPath(URI.create(s"oss://$bucket${fileName}/"))
    noException should be thrownBy ossProvider.checkAccess(dir)
  }

  it should "work for attribute view" taggedAs NeedAK in {
    val path = ossProvider.getPath(URI.create(s"oss://$bucket$fileName"))
    ossProvider.deleteIfExists(path)

    writeObject(path)
    val view = ossProvider.getFileAttributeView(path, classOf[OssStorageFileAttributesView])
    view shouldBe an [OssStorageFileAttributesView]

    val attr = view.readAttributes()
    attr shouldBe an [OssStorageObjectAttributes]

    val dir = ossProvider.getPath(URI.create(s"oss://$bucket${fileName}/"))
    val dirView = ossProvider.getFileAttributeView(dir, classOf[OssStorageFileAttributesView])
    dirView shouldBe an [OssStorageFileAttributesView]

    val dirAttr = dirView.readAttributes()
    dirAttr shouldBe an [OssStorageDirectoryAttributes]
  }

  it should "work for reading attrs" taggedAs NeedAK in {
    val path = ossProvider.getPath(URI.create(s"oss://$bucket$fileName"))
    ossProvider.deleteIfExists(path)

    writeObject(path)
    val attr = ossProvider.readAttributes(path, classOf[OssStorageFileAttributes])
    attr shouldBe an [OssStorageObjectAttributes]

    ossProvider.deleteIfExists(path)
    a [NoSuchFileException] should be thrownBy ossProvider.readAttributes(path, classOf[OssStorageFileAttributes])

    val dir = ossProvider.getPath(URI.create(s"oss://$bucket${fileName}/"))
    val dirAttr = ossProvider.readAttributes(dir, classOf[OssStorageFileAttributes])
    dirAttr shouldBe an [OssStorageDirectoryAttributes]
  }


  it should "work for reading dirs" taggedAs NeedAK in {
    val count = 10
    val testDir = "/test-read-dir"
    val filePrefix = "test-file"
    val expectedFileNames = ArrayBuffer.empty[String]
    val dir = ossProvider.getPath(URI.create(s"oss://$bucket$testDir/"))
    for (i <- 0 to count) {
      val fileName = filePrefix + i.toString
      expectedFileNames.append(fileName)

      val path = dir.resolve(fileName)

      ossProvider.deleteIfExists(path)
      writeObject(path)
    }

    val dirStream = ossProvider.newDirectoryStream(dir, new DirectoryStream.Filter[Path] {
      override def accept(entry: Path): Boolean = {
        true
      }
    })

    val files = ArrayBuffer.empty[String]
    dirStream.iterator.asScala foreach(file => files.append(file.toString))

    files should contain allElementsOf(expectedFileNames)
  }

}
