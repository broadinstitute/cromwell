package cromwell.engine.db

sealed trait CallBackendInfo

final case class LocalCallBackendInfo(processId: Option[Int]) extends CallBackendInfo

final case class JesCallBackendInfo(jesId: Option[JesId], jesStatus: Option[JesStatus]) extends CallBackendInfo

final case class SgeCallBackendInfo(sgeJobNumber: Option[Int]) extends CallBackendInfo
