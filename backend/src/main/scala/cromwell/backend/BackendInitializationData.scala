package cromwell.backend

// single-backend "marker interface"
trait BackendInitializationData

object BackendInitializationData {

  /**
    * Utility for extracting a type of backend data that we know is wrapped in an Option.
    * @param initializationDataOption The "optional" backend data passed usually from the initialization actor to each
    *                                 execution actor.
    * @tparam A The type to cast the initialization data.
    * @return The initialization data as the type A.
    */
  def as[A <: BackendInitializationData](initializationDataOption: Option[BackendInitializationData]): A = {
    initializationDataOption match {
      case Some(initializationData) => initializationData.asInstanceOf[A]
      case None => throw new RuntimeException("Initialization data was not found.")
    }
  }
}

object AllBackendInitializationData {
  def empty = AllBackendInitializationData(Map.empty)
}

// Holds initialization data for all backends initialized for a workflow.
case class AllBackendInitializationData(data: Map[String, Option[BackendInitializationData]]) {
  def get(backendName: String): Option[BackendInitializationData] = data.get(backendName).flatten
}
