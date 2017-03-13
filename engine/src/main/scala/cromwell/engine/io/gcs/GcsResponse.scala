package cromwell.engine.io.gcs

import cromwell.engine.io.IoActor._

import scala.language.existentials

/**
  * ADT used only inside the batch stream
  * @tparam T final type of the result of the Command
  */
private [gcs] sealed trait GcsBatchResponse[T]
private [gcs] case class GcsBatchTerminal[T](ioResult: IoResult) extends GcsBatchResponse[T]
private [gcs] case class GcsBatchRetry[T](context: GcsBatchCommandContext[T, _], failure: Throwable) extends GcsBatchResponse[T]
private [gcs] case class GcsBatchNextRequest[T](context: GcsBatchCommandContext[T, _]) extends GcsBatchResponse[T]
