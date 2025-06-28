package cromwell.backend.google.batch.util

import scala.util.{Failure, Success, Try}

sealed trait MemoryRetryCheckMode {
  val name: String

  override def toString: String = name
}

case object MemoryRetryRunnable extends MemoryRetryCheckMode {
  override val name: String = "Runnable"
}

case object MemoryRetryStandard extends MemoryRetryCheckMode {
  override val name: String = "Standard"
}

object MemoryRetryCheckMode {
  val DefaultMode = MemoryRetryRunnable
  private val AllModes: Seq[MemoryRetryCheckMode] = Seq(MemoryRetryRunnable, MemoryRetryStandard)

  def tryParse(mode: String): Try[MemoryRetryCheckMode] =
    AllModes
      .find(_.name.equalsIgnoreCase(mode))
      .map(Success(_))
      .getOrElse(
        Failure(
          new Exception(
            s"Invalid memory retry check mode: '$mode', supported modes are: ${AllModes.map(_.name).mkString("'", "', '", "'")}"
          )
        )
      )
}
