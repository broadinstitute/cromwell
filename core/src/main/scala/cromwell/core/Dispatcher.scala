package cromwell.core

object Dispatcher {
  val EngineDispatcher = "akka.dispatchers.engine-dispatcher"
  val IoDispatcher = "akka.dispatchers.io-dispatcher"
  val ApiDispatcher = "akka.dispatchers.api-dispatcher"
  val BackendDispatcher = "akka.dispatchers.backend-dispatcher"
  val ServiceDispatcher = "akka.dispatchers.service-dispatcher"
}
