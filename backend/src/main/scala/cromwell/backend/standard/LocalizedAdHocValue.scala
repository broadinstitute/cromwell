package cromwell.backend.standard

import cromwell.core.path.Path
import wom.callable.AdHocValue

/**
  * Represents an adhoc value that was moved to another location before the command is instantiated
  */
final case class LocalizedAdHocValue(originalValue: AdHocValue, localizedLocation: Path)
