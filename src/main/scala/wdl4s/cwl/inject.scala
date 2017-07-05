package shapeless
package syntax

object inject {
  import ops.coproduct.Inject

  /**
   * @author Fabio Labella
   */
  implicit class InjectSyntax[T](val t: T) extends AnyVal {
    /**
     * Inject the receiver into a coproduct `C`.
     * Only available if the coproduct contains the type `T`.
     */
    def inject[C <: Coproduct](implicit inj: Inject[C, T]): C = inj(t)
  }
}
