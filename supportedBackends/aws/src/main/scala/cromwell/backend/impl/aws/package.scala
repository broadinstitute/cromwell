package cromwell.backend.impl

import cats.data.ReaderT
import scala.language.higherKinds

package object aws {

  type Aws[F[_], A] = ReaderT[F, AwsBatchAttributes, A]

  def sanitize(name: String): String ={
    // Up to 128 letters (uppercase and lowercase), numbers, hyphens, and underscores are allowed.
    // We'll replace all invalid characters with an underscore
    name.replaceAll("[^A-Za-z0-9_\\-]", "_").slice(0, 128)
  }
}
