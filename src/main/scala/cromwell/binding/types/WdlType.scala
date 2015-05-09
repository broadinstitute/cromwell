package cromwell.binding.types

trait WdlType {
  def isCompatible(value: Any): Boolean

  def checkCompatible(value: Any) = if (!isCompatible(value)) throw new UnsupportedOperationException(s"Value '$value' is not compatible with WdlType $this.")
}
