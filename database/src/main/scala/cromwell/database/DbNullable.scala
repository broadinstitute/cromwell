package cromwell.database

/**
  * Trait and typeclasses for converting to and from non-null representations of data in the database that is
  * actually null in the domain model.
  * This is used to work around https://bugs.mysql.com/bug.php?id=8173
  */
object DbNullable {

  trait DbNullable[A] {
    def toDb(in: Option[A]): A
    def fromDb(out: A): Option[A]
  }

  implicit object DbNullableInt extends DbNullable[Int] {
    val NoInt = -1
    override def toDb(in: Option[Int]): Int = in getOrElse NoInt
    override def fromDb(out: Int): Option[Int] = if (out == NoInt) None else Option(out)
  }

  implicit object DbNullableString extends DbNullable[String] {
    val NoString = ""
    override def toDb(in: Option[String]): String = in getOrElse NoString
    override def fromDb(out: String): Option[String] = if (out == NoString) None else Option(out)
  }

  implicit class ToDbNullableEnhancer[A](val a: Option[A])(implicit N: DbNullable[A]) {
    def toDb: A = N.toDb(a)
  }

  implicit class FromDbNullableEnhancer[A](val a: A)(implicit N: DbNullable[A]) {
    def fromDb: Option[A] = N.fromDb(a)
  }
}
