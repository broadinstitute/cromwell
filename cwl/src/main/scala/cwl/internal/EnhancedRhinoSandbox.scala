package cwl.internal

import cwl.internal.EnhancedRhinoSandbox._
import delight.rhinosandox.internal._
import org.mozilla.javascript._

import scala.jdk.CollectionConverters._
import scala.reflect._

/**
  * Extends the RhinoSandboxImpl with some fixes and enhancements using java reflection.
  *
  * @param strict Should evaluation be strict.
  * @param languageVersionOption The optional language version to set.
  */
class EnhancedRhinoSandbox(strict: Boolean = true, languageVersionOption: Option[Int] = None) extends RhinoSandboxImpl {

  // Allows easier reflection access to private fields
  private lazy val sandboxImpl: RhinoSandboxImpl = this

  /**
    * Similar to RhinoSandbox.eval but passes back the context and scope for mutating before evaluation.
    *
    * With the original RhinoSandbox.eval not sure how to:
    * - Create (nested) arrays
    * - Create (nested) maps
    * - Set JS/ES version
    *
    * So this uses a copy-port of the original, using reflection to read some of RhinoSandbox's private variables.
    *
    * Instead of hiding the context and scope as in RhinoSandbox.eval, both are passed back through the block.
    *
    * TODO: Ask RhinoSandbox if hacks are even needed, and if so contrib back patches so that reflection isn't required.
    * - Is there a way to pass nested maps/arrays via a `java.util.Map[String, Object]`, or must we use our `block`?
    * - Additionally: can we skip passing `context` to a block as a thread may just call `Context.getCurrentContext()`?
    *
    * @see https://maxrohde.com/2015/08/06/sandboxing-javascript-in-java-app-link-collection/
    * @see https://github.com/javadelight/delight-rhino-sandbox/blob/9f5a073/src/main/java/delight/rhinosandox/internal/RhinoSandboxImpl.java#L100-L123
    * @see delight.rhinosandox.internal.RhinoSandboxImpl#assertContextFactory()
    */
  def eval(sourceName: String, js: String)(block: (Context, Scriptable) => Unit): AnyRef = {
    assertContextFactory()
    val sandboxImpl_contextFactory = PrivateField(sandboxImpl, "contextFactory").as[ContextFactory]
    // RhinoSandbox diff: eval has enterContext inside the try, but Rhino docs say it belongs outside.
    // https://www-archive.mozilla.org/rhino/apidocs/org/mozilla/javascript/ContextFactory.html#enterContext%28%29
    val context = sandboxImpl_contextFactory.enterContext
    try {
      // RhinoSandbox diff: allow setting the language version
      languageVersionOption foreach context.setLanguageVersion
      assertSafeScope(context)
      val sandboxImpl_globalScope = PrivateField(sandboxImpl, "globalScope").as[ScriptableObject]
      val sandboxImpl_sealScope = PrivateField(sandboxImpl, "sealScope").as[Boolean]
      if (sandboxImpl_sealScope) {
        sandboxImpl_globalScope.sealObject()
      }
      val sandboxImpl_safeScope = PrivateField(sandboxImpl, "safeScope").as[ScriptableObject]
      val instanceScope = context.newObject(sandboxImpl_safeScope)
      instanceScope.setPrototype(sandboxImpl_safeScope)
      instanceScope.setParentScope(null)

      block(context, instanceScope)

      // RhinoSandbox diff: allow strict JS/ES evaluation
      // See note at top assertContextFactory as to why we have to put 'use strict'; here.
      // All on one line to avoid off-by-one error for javascript error messages that report line numbers.
      // Could also pass zero as the line number, but the RhinoSandbox passes hard codes line number one also.
      val script = if (strict) s"'use strict';$js" else js
      context.evaluateString(instanceScope, script, sourceName, 1, null)
    } finally {
      Context.exit()
    }
  }

  /**
    * Stricter context factory modified from RhinoSandboxImpl.assertContextFactory().
    *
    * ContextFactory.initGlobal() is called within a static synchronized block.
    * The globalScope initialized via initSafeStandardObjects instead of initStandardObjects.
    *
    * The default implementation uses a SafeContext that allows non-strict JS/ES. We would ideally set
    * FEATURE_STRICT_MODE to true but that only produces warnings and doesn't return an error. Unfortunately when
    * FEATURE_WARNING_AS_ERROR is enabled then non-strict Rhino warnings like "missing ;" throw errors. Instead,
    * "'use strict';" is injected before scripts.
    */
  override def assertContextFactory(): Unit = {
    if (PrivateField(sandboxImpl, "contextFactory").as[ContextFactory] != null) {
      return
    }

    val _safeContext = new SafeContext
    PrivateField(sandboxImpl, "contextFactory") := _safeContext
    val _hasExplicitGlobal = ContextFactory.hasExplicitGlobal
    val _not = !_hasExplicitGlobal
    // RhinoSandbox diff: the global does not like to be initialized twice. Synchronize initialization.
    val sandboxImpl_contextFactory = PrivateField(sandboxImpl, "contextFactory").as[SafeContext]
    if (_not) initGlobalSynchronized(sandboxImpl_contextFactory)

    val sandboxImpl_instructionLimit = PrivateField(sandboxImpl, "instructionLimit").as[Int]
    PrivateField(sandboxImpl_contextFactory, "maxInstructions") := sandboxImpl_instructionLimit
    val sandboxImpl_maxDuration = PrivateField(sandboxImpl, "maxDuration").as[Long]
    PrivateField(sandboxImpl_contextFactory, "maxRuntimeInMs") := sandboxImpl_maxDuration
    // RhinoSandbox diff: assertContextFactory has enterContext inside the try, but Rhino docs say it belongs outside.
    // https://www-archive.mozilla.org/rhino/apidocs/org/mozilla/javascript/ContextFactory.html#enterContext%28%29
    val context = sandboxImpl_contextFactory.enterContext
    try {
      // RhinoSandbox diff: Default globalScope is created via initStandardObjects instead of initSafeStandardObjects.
      // initStandardObjects would add the various java packages into the global scope, including `java.io.File`, etc.
      PrivateField(sandboxImpl, "globalScope") := context.initSafeStandardObjects(null, false)
      val sandboxImpl_inScope = PrivateField(sandboxImpl, "inScope").as[java.util.Map[String, AnyRef]]
      val _entrySet = sandboxImpl_inScope.entrySet
      val sandboxImpl_globalScope = PrivateField(sandboxImpl, "globalScope").as[ScriptableObject]
      for (entry <- _entrySet.asScala) {
        sandboxImpl_globalScope.put(
          entry.getKey,
          sandboxImpl_globalScope,
          Context.toObject(entry.getValue, sandboxImpl_globalScope))
      }
      val parameters = Array(classOf[String])
      val dealMethod = classOf[RhinoEvalDummy].getMethod("eval", parameters: _*)
      val _rhinoEval = new RhinoEval("eval", dealMethod, sandboxImpl_globalScope)
      sandboxImpl_globalScope.defineProperty("eval", _rhinoEval, ScriptableObject.DONTENUM)
    } finally {
      Context.exit()
    }
  }

}

object EnhancedRhinoSandbox {

  /**
    * Get or sets a private field.
    *
    * @param obj  The instance to retrieve the field from.
    * @param name The name of the field.
    * @tparam A The class to retrieve the value from. The field MUST exist on this class, and not a superclass.
    */
  final case class PrivateField[A: ClassTag](obj: A, name: String) {
    private[this] def field = {
      val field = classTag[A].runtimeClass.getDeclaredField(name)
      field.setAccessible(true)
      field
    }

    def as[B]: B = {
      field.get(obj).asInstanceOf[B]
    }

    def :=(value: Any): Unit = {
      field.set(obj, value)
    }
  }

  /**
    * Call ContextFactory.initGlobal in a static synchronized block.
    *
    * @see [[cwl.internal.EnhancedRhinoSandbox.assertContextFactory]]
    */
  private def initGlobalSynchronized(sandboxImpl_contextFactory: ContextFactory) = {
    synchronized {
      val _hasExplicitGlobal = ContextFactory.hasExplicitGlobal
      val _not = !_hasExplicitGlobal
      if (_not) ContextFactory.initGlobal(sandboxImpl_contextFactory)
    }
  }

}
