package cromwell.engine.workflow

package object lifecycle {
  case object EngineLifecycleActorAbortCommand

  trait EngineLifecycleStateCompleteResponse
  trait EngineLifecycleActorAbortedResponse extends EngineLifecycleStateCompleteResponse
}
