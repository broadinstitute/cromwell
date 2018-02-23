package wdl.draft3.transforms.wdlom2wom.linking

sealed trait UnlinkedGeneratedValueName { def linkableName: String }

sealed trait UnlinkedConsumedValueName

final case class UnlinkedIdentifierName(name: String) extends UnlinkedConsumedValueName with UnlinkedGeneratedValueName {
  override def linkableName = name
}

final case class UnlinkedCallOutputName(taskName: String, outputName: String) extends UnlinkedGeneratedValueName {
  override def linkableName: String = s"$taskName.$outputName"
}

/**
  * Until we do the linking, we can't tell whether a consumed 'x.y' is a call output or a member access for 'y' on
  * a variable called 'x'.
  */
final case class UnlinkedCallOutputOrIdentifierAndMemberAccess(name: String, firstLookup: String) extends UnlinkedConsumedValueName
