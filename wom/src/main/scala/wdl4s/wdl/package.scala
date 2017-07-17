package wdl4s.wdl

import org.apache.commons.codec.digest.DigestUtils

import scala.util.Try

trait TsvSerializable {
  def tsvSerialize: Try[String]
}

package object values {

  type FileHasher = WdlFile => SymbolHash

  implicit class HashableString(val value: String) extends AnyVal with Hashable {
    def md5Sum: String = DigestUtils.md5Hex(value)
    def md5SumShort: String = value.md5Sum.substring(0, 8)
  }
}
