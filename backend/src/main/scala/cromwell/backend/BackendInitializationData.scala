package cromwell.backend

// single-backend "marker interface"
trait BackendInitializationData

object AllBackendInitializationData {
  def empty = AllBackendInitializationData(Map.empty)
}

// Holds initialization data for all backends initialized for a workflow.
case class AllBackendInitializationData(data: Map[String, Option[BackendInitializationData]]) {
  def get(backendName: String): Option[BackendInitializationData] = data.get(backendName).flatten
}
