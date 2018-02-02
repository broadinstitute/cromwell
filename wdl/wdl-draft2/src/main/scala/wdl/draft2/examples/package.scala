package wdl.draft2

import wdl.versioning.WdlVersionSpecifics

package object examples {
  implicit val wdlVersionSpecifics: WdlVersionSpecifics = Draft2VersionSpecifics
}
