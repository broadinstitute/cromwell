package centaur.test.workflow

sealed trait BackendsRequirement

object BackendsRequirement {
  def fromConfig(backendMode: String, backendsList: List[String]): BackendsRequirement = backendMode match {
    case "all" => AllBackendsRequired(backendsList)
    case "only" => OnlyBackendsAllowed(backendsList)
    case _ => AnyBackendRequired(backendsList)
  }
}

final case class AllBackendsRequired(backends: List[String]) extends BackendsRequirement
final case class AnyBackendRequired(backends: List[String]) extends BackendsRequirement
final case class OnlyBackendsAllowed(backends: List[String]) extends BackendsRequirement
