package cromwell.core.instrumentation

import cats.data.NonEmptyList

object InstrumentationKeys {
  val SuccessKey = NonEmptyList.of("success")
  val AbortedKey = NonEmptyList.of("aborted")
  val FailureKey = NonEmptyList.of("failure")
  val RetryKey = NonEmptyList.of("retry")
  val CountKey = NonEmptyList.of("count")
  val TimingKey = NonEmptyList.of("timing")
}
