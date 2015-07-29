package cromwell.binding

// FIXME: It'd be nice to have a notion of a parented and a not-parented Scope, see FIXME in Call about this
trait Scope {
  def name: String

  def parent: Option[Scope]

  def fullyQualifiedName: String = parent match {
    case Some(x: Scope) => s"${x.fullyQualifiedName}.$name"
    case _ => name
  }

}
