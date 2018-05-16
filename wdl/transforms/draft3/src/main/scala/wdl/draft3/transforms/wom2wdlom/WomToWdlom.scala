package wdl.draft3.transforms.wom2wdlom

trait WomToWdlom[A, B] {
  def convert(a: A): B
}

object WomToWdlom {
  // https://stackoverflow.com/a/48915864/818054
  def apply[A, B](implicit converter: WomToWdlom[A, B]): WomToWdlom[A, B] = converter

  object ops {
    def convert[A, B](a: A)(implicit converter: WomToWdlom[A, B]): B = converter.convert(a)

    implicit class WomToWdlomOps[A, B](a: A)(implicit converter: WomToWdlom[A, B]) {
      def convert(implicit converter: WomToWdlom[A, B]): B = converter.convert(a)
    }
  }
}
