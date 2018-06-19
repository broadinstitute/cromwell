package womtool.wom2wdlom

trait WomToWdlom[A, B] {
  def toWdlom(a: A): B
}

object WomToWdlom {
  // https://stackoverflow.com/a/48915864/818054
  def apply[A, B](implicit converter: WomToWdlom[A, B]): WomToWdlom[A, B] = converter

  object ops {
    def toWdlom[A, B](a: A)(implicit converter: WomToWdlom[A, B]): B = converter.toWdlom(a)

    implicit class WomToWdlomOps[A, B](a: A)(implicit converter: WomToWdlom[A, B]) {
      def toWdlom(implicit converter: WomToWdlom[A, B]): B = converter.toWdlom(a)
    }
  }
}
