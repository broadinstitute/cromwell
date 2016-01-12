package cromwell.backend.model

/**
  * Defines event types.
  */
trait EventType
/**
  * Represents an event due to a task execution.
  */
case object ExecutionEvent extends EventType

/**
  * Serves as a subscription object in which the subscriber is tied to an event.
  * @param eventType Event type.
  * @param subscriber Subscriber to that kind of event. It can be an Akka actor or any other sort of implementation.
  * @tparam T Subscriber type.
  */
case class Subscription[T](eventType: EventType, subscriber: T)