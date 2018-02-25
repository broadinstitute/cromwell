package centaur.api

import java.util.concurrent.ThreadFactory
import java.util.concurrent.atomic.AtomicInteger

/**
  * A static version of java.util.concurrent.Executors.DefaultThreadFactory that creates daemon threads that exit when
  * the application exits.
  *
  * See also: java.util.concurrent.Executors.DefaultThreadFactory
  * See also: http://dev.bizo.com/2014/06/cached-thread-pool-considered-harmlful.html
  * > Always provide your own ThreadFactory so all your threads are named appropriately and have the daemon flag set.
  * > Your thread pool shouldn't keep the application alive. That's the responsibility of the application itself (i.e.,
  * > main thread).
  */
object DaemonizedDefaultThreadFactory extends ThreadFactory {
  private val s = System.getSecurityManager
  private val group = if (s != null) s.getThreadGroup else Thread.currentThread.getThreadGroup
  private val threadNumber = new AtomicInteger(1)
  private val namePrefix = "daemonpool-thread-"

  override def newThread(r: Runnable): Thread = {
    val t = new Thread(group, r, namePrefix + threadNumber.getAndIncrement, 0)
    if (!t.isDaemon) t.setDaemon(true)
    if (t.getPriority != Thread.NORM_PRIORITY) t.setPriority(Thread.NORM_PRIORITY)
    t
  }
}
