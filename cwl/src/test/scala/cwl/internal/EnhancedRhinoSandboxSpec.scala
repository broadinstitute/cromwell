package cwl.internal

import java.util.concurrent.{CountDownLatch, TimeUnit}

import cwl.internal.EnhancedRhinoSandboxSpec.CachedClassLoader
import org.apache.commons.io.IOUtils
import org.mozilla.javascript.ContextFactory
import org.scalatest.concurrent.ScaledTimeSpans
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.concurrent.duration._

class EnhancedRhinoSandboxSpec extends AnyFlatSpec with Matchers with ScaledTimeSpans {

  behavior of "EnhancedRhinoSandbox"

  it should "synchronize global ContextFactory initialization" in {
    // Initialize the context factory here to be sure that our loaded runner is actually using a fresh ContextFactory
    new EnhancedRhinoSandbox(true).assertContextFactory()

    val threadCount = 10

    // A custom loader that will contain our fresh static ContextFactory
    val cachedClassLoader = new CachedClassLoader
    // A latch that lets this test know the threads are ready
    val threadReadyLatch = new CountDownLatch(threadCount)
    // A latch to ensure that all threads try to run at the same time
    val threadGoLatch = new CountDownLatch(1)

    class Runner(name: String) extends Thread {
      var success = false

      override def run(): Unit = {
        EnhancedRhinoSandboxSpec.runClassLoaded(
          cachedClassLoader,
          threadReadyLatch = threadReadyLatch,
          threadGoLatch = threadGoLatch
        )
        success = true
      }
    }

    val threads = (0 until threadCount).map(i => new Runner(i.toString))
    // move all the threads to the starting line...
    threads.foreach(_.start())
    // Wait for threads to say they're ready
    threadReadyLatch.await(scaled(10.seconds).toMillis, TimeUnit.MILLISECONDS)
    // go!
    threadGoLatch.countDown()
    // Maximum time to wait for each runner to finish
    threads.foreach(_.join(scaled(10.seconds).toMillis))

    withClue("all threads did not complete successfully:") {
      threads.map(_.success) should contain only true
    }
  }
}

object EnhancedRhinoSandboxSpec extends ScaledTimeSpans {
  /**
    * Runs a runner loaded from a custom class loader.
    */
  private def runClassLoaded(ecmaScriptUtilClassLoader: CachedClassLoader,
                             threadReadyLatch: CountDownLatch,
                             threadGoLatch: CountDownLatch): Unit = {
    val cls = ecmaScriptUtilClassLoader.loadClass(classOf[LatchedRunner].getName)
    val obj = cls.newInstance
    cls
      .getMethod("run", classOf[CountDownLatch], classOf[CountDownLatch])
      .invoke(obj, threadReadyLatch, threadGoLatch)
    ()
  }

  /**
    * Waits for a shared latch and then quickly tries to run the EcmaScriptUtil.
    * Before running, there is a test that ensures that the ContextFactory has _not_ been initialized.
    *
    * Inspired by: https://stackoverflow.com/questions/13263079/run-two-thread-at-the-same-time-in-java#13263273
    */
  class LatchedRunner {
    def run(threadReadyLatch: CountDownLatch, threadGoLatch: CountDownLatch): Unit = {
      require(!ContextFactory.hasExplicitGlobal, "The ContextFactory within this classloader is already initialized")
      // Notify the test we're at the starting line.
      threadReadyLatch.countDown()
      // This await is for the time between all the threads `start()`-ing up, and the `countDownLatch.countDown()`
      threadGoLatch.await(scaled(5.seconds).toMillis, TimeUnit.MILLISECONDS)
      new EnhancedRhinoSandbox(true).assertContextFactory()
      ()
    }
  }

  /**
    * Locates a class but instead of returning the class from the original classloader it instead returns a cached copy.
    *
    * This is needed when we want to ensure that a static variable has NOT been initialized on a class instance.
    *
    * Parent class loader must be the current thread's class loader since the system class loader is NOT sufficient
    * under sbt.
    *
    * Inspired by: https://stackoverflow.com/questions/14316996/how-to-replace-classes-in-a-running-application-in-java#14317506
    */
  class CachedClassLoader extends ClassLoader(Thread.currentThread.getContextClassLoader) {

    override def loadClass(name: String): Class[_] = {
      loadedClasses synchronized {
        loadedClasses.getOrElseUpdate(name, internalLoadClass(name))
      }
    }

    /** Classes to always load from the parent. */
    private val parentPackages = List("scala", "java", "sun")
    /** Synchronized mutable map of our already loaded classes. */
    private val loadedClasses = scala.collection.mutable.Map.empty[String, Class[_]]

    /**
      * Locates a class but returns our _own_ copy of it. Classes under "scala", "java", and "sun" will use the parent.
      * If a class cannot be located for any reason, the parent is also returned.
      *
      * @param name Name of the class to load.
      * @return The loaded class.
      */
    private def internalLoadClass(name: String): Class[_] = {
      if (!parentPackages.exists(pkg => name.startsWith(pkg + "."))) {
        val inputStream = getParent.getResourceAsStream(name.replace(".", "/").concat(".class"))
        try {
          val bytes = IOUtils.toByteArray(inputStream)
          return defineClass(name, bytes, 0, bytes.length)
        } catch {
          case _: Throwable => /* ignore and fall through to loading from the parent */
        } finally {
          if (inputStream != null)
            inputStream.close()
        }
      }
      getParent.loadClass(name)
    }
  }

}
