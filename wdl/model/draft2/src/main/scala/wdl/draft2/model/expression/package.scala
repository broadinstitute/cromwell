package wdl.draft2.model

import wdl.shared.FileSizeLimitationConfig

import scala.util.Try

package object expression {

  type WdlStandardLibraryReadFunction[A] = FileSizeLimitationConfig => Try[A]

}
