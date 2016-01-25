package cromwell.backend.model

/**
  * Defines event types.
  */
trait EventType

/**
  * Represents an event due to a task execution.
  */
class ExecutionEvent extends EventType

trait SubscriptionEvent extends EventType

case object Subscribed extends SubscriptionEvent

case object Unsubscribed extends SubscriptionEvent

/**
  * Serves as a subscription object in which the subscriber is tied to an event.
  * @param eventType Event type.
  * @param subscriber Subscriber to that kind of event. It can be an Akka actor or any other sort of implementation.
  * @tparam A Subscriber type.
  */
case class Subscription[A](eventType: EventType, subscriber: A)