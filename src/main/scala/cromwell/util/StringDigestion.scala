package cromwell.util

import java.io.{File, FileInputStream}

import org.apache.commons.codec.digest.DigestUtils

object StringDigestion {
  implicit class StringDigester(val string: String) extends AnyVal {
    def md5Sum: String = DigestUtils.md5Hex(string)
  }

  implicit class FileDigester(val file: File) extends AnyVal {
    def md5Sum: String = {
      val fis = new FileInputStream(file)
      try {
        DigestUtils.md5Hex(fis)
      } finally fis.close()
    }
  }
}