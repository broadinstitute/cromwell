package wdl4s

import org.apache.commons.codec.digest.DigestUtils


package object values {

  case class SymbolHash(value: String) extends AnyVal with Ordered[SymbolHash] {
    def compare(that: SymbolHash) = this.value compare that.value
  }

  trait Hashable extends Any {
    def md5Sum: String
  }

  type FileHasher = WdlFile => SymbolHash

  implicit class HashableString(val value: String) extends AnyVal with Hashable {
    def md5Sum: String = DigestUtils.md5Hex(value)
    def md5SumShort: String = value.md5Sum.substring(0, 8)
  }
}
