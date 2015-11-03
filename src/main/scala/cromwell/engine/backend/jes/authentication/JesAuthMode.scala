package cromwell.engine.backend.jes.authentication

object JesAuthMode {
  def fromString(name: String): JesAuthMode = name match {
    case "service" => ServiceAccountMode
    case "user" => UserMode
    case nop => throw new IllegalArgumentException(s"$nop is not a recognized authentication mode")
  }
}

// Authentication modes supported
sealed trait JesAuthMode
object ServiceAccountMode extends JesAuthMode
object UserMode extends JesAuthMode
