package cromwell.backend.impl

import cats.data.ReaderT
import scala.language.higherKinds

package object aws {

  type Aws[F[_], A] = ReaderT[F, AwsBatchAttributes, A]
}
