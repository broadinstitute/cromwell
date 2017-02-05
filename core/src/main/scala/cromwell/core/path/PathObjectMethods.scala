package cromwell.core.path

/**
  * Overrides java.lang.Object methods for Path
  */
trait PathObjectMethods {
  self: Path =>

  override def toString: String = pathAsString

  override def equals(obj: Any) = {
    obj match {
      case other: Path => nioPathPrivate == other.nioPathPrivate
      case _ => false
    }
  }

  override def hashCode = nioPathPrivate.hashCode()
}
