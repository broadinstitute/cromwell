package cloud.nio.util

import cats.data.NonEmptyList

class AggregateException(message: String, val exceptions: NonEmptyList[Exception]) extends Exception(message) {
  exceptions.toList.foreach(addSuppressed)
}
