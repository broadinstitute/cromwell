package cromwell.backend

import cromwell.core.simpleton.WomValueSimpleton

object BackendCacheHitCopyingActor {
  final case class CopyOutputsCommand(womValueSimpletons: Seq[WomValueSimpleton], jobDetritusFiles: Map[String, String], returnCode: Option[Int])
}
