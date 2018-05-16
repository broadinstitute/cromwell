package wdl.draft3.transforms.wom2wdlom

trait Convert[A, B] {
  def convert(a: A): B
}

object Convert {
  // https://stackoverflow.com/a/48915864/818054
  def apply[A, B](implicit converter: Convert[A, B]): Convert[A, B] = converter

  object ops {
    def convert[A, B](a: A)(implicit converter: Convert[A, B]): B = converter.convert(a)

    implicit class ConvertOps[A, B](a: A)(implicit converter: Convert[A, B]) {
      def convert(implicit converter: Convert[A, B]): B = converter.convert(a)
    }
  }
}
