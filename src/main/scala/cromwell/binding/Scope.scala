package cromwell.binding

trait Scope {
  def name: String

  def parent: Option[Scope]

  def fullyQualifiedName: String = parent match {
    case Some(x: Scope) => s"${x.fullyQualifiedName}.$name"
    case _ => name
  }

}
