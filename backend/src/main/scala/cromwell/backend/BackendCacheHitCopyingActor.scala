package cromwell.backend

import cromwell.core.simpleton.WdlValueSimpleton

object BackendCacheHitCopyingActor {
  final case class CopyOutputsCommand(wdlValueSimpletons: Seq[WdlValueSimpleton], jobDetritusFiles: Map[String, String], returnCode: Option[Int])
}
