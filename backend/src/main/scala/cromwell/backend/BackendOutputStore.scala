package cromwell.backend

import wdl4s.Scope
import wdl4s.types.WdlType
import wdl4s.values.WdlValue

// TODO: PBE: Not moving the "engine" version yet as it has many other "engine" classes mixed in. Eventually "core"?

case class BackendOutputEntry(name: String, wdlType: WdlType, wdlValue: Option[WdlValue])

case class BackendOutputKey(scope: Scope, index: Option[Int])

case class BackendOutputStore(store: Map[BackendOutputKey, Traversable[BackendOutputEntry]])
