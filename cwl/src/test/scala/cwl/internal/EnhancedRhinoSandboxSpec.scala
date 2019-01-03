package cwl.internal

import java.util.concurrent.{CountDownLatch, TimeUnit}

import cwl.internal.EnhancedRhinoSandboxSpec.CachedClassLoader
import org.apache.commons.io.IOUtils
import org.mozilla.javascript.ContextFactory
import org.scalatest.{FlatSpec, Matchers}

class EnhancedRhinoSandboxSpec extends FlatSpec with Matchers {

  behavior of "EnhancedRhinoSandbox"

  it should "synchronize global ContextFactory initialization" in {
    // Initialize the context factory here to be sure that our loaded runner is actually using a fresh ContextFactory
    new EnhancedRhinoSandbox(true).assertContextFactory()

    // A custom loader that will contain our fresh static ContextFactory
    val cachedClassLoader = new CachedClassLoader
    // A latch to ensure that all threads try to run at the same time
    val countDownLatch = new CountDownLatch(1)

    class Runner(name: String) extends Thread {
      var success = false

      override def run(): Unit = {
        EnhancedRhinoSandboxSpec.runClassLoaded(cachedClassLoader, countDownLatch)
        success = true
      }
    }

    val threads = (0 until 10).map(i => new Runner(i.toString))
    threads.foreach(_.start())
    countDownLatch.countDown()
    threads.foreach(_.join(1000L))
    withClue("all threads did not complete successfully:") {
      threads.map(_.success) should contain only true
    }
  }
}

object EnhancedRhinoSandboxSpec {
  /**
    * Runs a runner loaded from a custom class loader.
    */
  private def runClassLoaded(ecmaScriptUtilClassLoader: CachedClassLoader,
                             countDownLatch: CountDownLatch): Unit = {
    val cls = ecmaScriptUtilClassLoader.loadClass(classOf[LatchedRunner].getName)
    val obj = cls.newInstance
    cls.getMethod("run", classOf[CountDownLatch]).invoke(obj, countDownLatch)
    ()
  }

  /**
    * Waits for a shared latch and then quickly tries to run the EcmaScriptUtil.
    * Before running, there is a test that ensures that the ContextFactory has _not_ been initialized.
    *
    * Inspired by: https://stackoverflow.com/questions/13263079/run-two-thread-at-the-same-time-in-java#13263273
    */
  class LatchedRunner {
    def run(countDownLatch: CountDownLatch): Unit = {
      require(!ContextFactory.hasExplicitGlobal, "The ContextFactory within this classloader is already initialized")
      countDownLatch.await(5, TimeUnit.SECONDS)
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
