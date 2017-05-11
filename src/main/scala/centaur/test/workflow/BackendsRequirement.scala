package centaur.test.workflow

sealed trait BackendsRequirement

object BackendsRequirement {
  def fromConfig(backendMode: String, backendsList: List[String]): BackendsRequirement = if (backendMode == "all") {
    AllBackendsRequired(backendsList)
  } else {
    AnyBackendRequired(backendsList)
  }
}

final case class AllBackendsRequired(backends: List[String]) extends BackendsRequirement
final case class AnyBackendRequired(backends: List[String]) extends BackendsRequirement
