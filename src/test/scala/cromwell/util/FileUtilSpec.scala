package cromwell.util

import org.scalatest.{Matchers, FlatSpec}
import cromwell.util.FileUtil._

class FileUtilSpec extends FlatSpec with Matchers with TestFileUtil {

  it should "generate a md5 hash from a java.io.File" in {
    val f1 = createCannedFile("file1", "content")
    val f2 = createCannedFile("file2", "content")
    val f3 = createCannedFile("file3", "other content")
    f1.md5Sum should be(f2.md5Sum)
    f1.md5Sum shouldNot be(f3.md5Sum)
  }

  it should "slurp content from a java.io.File" in {
    val f1 = createCannedFile("file1", "content")
    f1.slurp should be("content")
  }

}
