package cromwell.util

import org.scalatest.{Matchers, FlatSpec}
import cromwell.util.FileUtil._

class FileUtilSpec extends FlatSpec with Matchers with TestFileUtil {

  it should "generate a md5 hash from a file" in {
    val f1 = createCannedFile("file1", "content")
    val f2 = createCannedFile("file2", "content")
    val f3 = createCannedFile("file3", "other content")
    f1.md5Sum should be("9a0364b9e99bb480dd25e1f0284c8555")
    f1.md5Sum should be(f2.md5Sum)
    f1.md5Sum shouldNot be(f3.md5Sum)
  }
}
