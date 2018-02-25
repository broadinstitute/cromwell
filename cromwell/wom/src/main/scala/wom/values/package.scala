package wom.values

case class SymbolHash(value: String) extends AnyVal with Ordered[SymbolHash] {
  def compare(that: SymbolHash) = this.value compare that.value
}

object GlobFunctions {
  def globName(glob: String) = s"glob-${glob.md5Sum}"
  def prefixWithGlobDir(glob: String) = s"${globName(glob)}/$glob"
}

trait Hashable extends Any {
  def md5Sum: String
}
