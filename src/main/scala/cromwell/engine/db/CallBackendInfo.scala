package cromwell.engine.db

sealed trait CallBackendInfo

final case class LocalCallBackendInfo(processId: Option[Int]) extends CallBackendInfo
