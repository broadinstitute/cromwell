package cromwell.backend.impl.htcondor.caching

import akka.actor.Props

trait CacheActorFactory {

  def getCacheActorProps(forceRewrite: Boolean): Props

}
