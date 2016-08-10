package cromwell.util

import scala.util.Try

/**
  * If this looks java-y, it's because it is. It's (almost) a direct-to-scala mapping of the Java try-with-resources block.
  * This code looks ugly so yours won't have to! (TM)
  */
object TryWithResource {
  def tryWithResource[A <: AutoCloseable, B](getResource: () => A)(block: A => B) = Try {
    var t: Option[Throwable] = None
    var resource: Option[A] = None
    try {
      resource = Option(getResource())
      resource match {
        case Some(r) => block(r)
        case None => throw new NullPointerException("Expected resource could not be created")
      }
    } catch {
      case x: Throwable =>
        t = Option(x)
        throw x
    } finally {
      resource foreach { r =>
        try {
          r.close()
        } catch {
          case y: Throwable => t match {
            case Some(_t) => _t.addSuppressed(y)
            case None => throw y
          }
        }
      }
    }
  }
}
