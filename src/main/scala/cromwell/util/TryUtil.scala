package cromwell.util

import java.io.{PrintWriter, StringWriter}

import scala.util.{Failure, Try}

object TryUtil {

  def stringifyFailures[T](failure: Try[T]): String = {
    val writer = new PrintWriter(new StringWriter())
    failure.recover { case e => e.printStackTrace(writer)}
    writer.flush()
    writer.close()
    writer.toString
  }

  def stringifyFailures[T](possibleFailures: Traversable[Try[T]]): Traversable[String] =
    possibleFailures.collect { case failure: Failure[T] => stringifyFailures(failure) }
}
