package cromwell.languages.util


class ResolvedImportsStore {
  private val mutex = new Object
  private var resolvedImportsSet = Set.empty[String]

  def updateResolvedImportsSet(importPath: String): Unit = {
    mutex.synchronized {
      resolvedImportsSet = resolvedImportsSet + importPath
    }
  }

  def getResolvedImportsSet: Set[String] = resolvedImportsSet
}


object NoopResolvedImportsStore extends ResolvedImportsStore {
  override def updateResolvedImportsSet(importPath: String): Unit = {}
}
