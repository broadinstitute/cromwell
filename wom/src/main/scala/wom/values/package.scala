package wom.values

case class SymbolHash(value: String) extends AnyVal with Ordered[SymbolHash] {
  def compare(that: SymbolHash) = this.value compare that.value
}

trait Hashable extends Any {
  def md5Sum: String
}
