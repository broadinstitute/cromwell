package common.mock

import org.mockito.{ArgumentCaptor, Mockito}

import scala.reflect.{classTag, ClassTag}

/**
  * Yet another scala wrapper around Mockito.
  *
  * As of Aug 2021 there are a few mockito wrapper choices.
  *
  * `mockito-scala`:
  *  - lots of nice DSL, but who knows if it'll last
  *  - behind on Scala 3 support
  *    - [[https://github.com/mockito/mockito-scala/issues/364]]
  *
  * `scalatestplus`:
  *  - stuck on mockito 3.4.x
  *    - [[https://github.com/scalatest/scalatestplus-mockito/issues/24]]
  *  - entire library only provides four one-line wrappers
  *    - [[https://github.com/scalatest/scalatestplus-mockito/blob/release-3.2.9.0-for-mockito-3.4/src/main/scala/org/scalatestplus/mockito/MockitoSugar.scala]]
  *
  * `scalamock`:
  *  - might be abandoned?
  *    - [[https://github.com/paulbutcher/ScalaMock/issues/396]]
  *
  * `specs2-mock`:
  *  - As of Specs2 5.x the mock wrappers appear to be gone, pointing to scalamock instead:
  *    - [[https://etorreborre.github.io/specs2/guide/5.0.0-RC-01/org.specs2.guide.Installation.html]]
  *  - Btw, specs2-mock 4.x pulls in specs2-core possibly leading SBT to look for specs using the specs2 test framework:
  *    - [[https://etorreborre.github.io/specs2/guide/SPECS2-4.12.0/org.specs2.guide.Installation.html#other-dependencies]]
  */
trait MockSugar extends MockImplicits {

  /**
    * Returns a new mock with Smart Nulls.
    *
    * Note: if you run into issues with `mock` then try [[mockWithDefaults]].
    */
  def mock[A: ClassTag]: A =
    Mockito.mock(
      classTag[A].runtimeClass.asInstanceOf[Class[A]],
      Mockito.withSettings().defaultAnswer(Mockito.RETURNS_SMART_NULLS)
    )

  /**
    * Creates a mock returning default values instead of Smart Nulls.
    *
    * Works around a cryptic issue that that popped up in PipelinesApiBackendCacheHitCopyingActorSpec:
    * {{{
    * Underlying exception : java.lang.IllegalArgumentException: Cannot cast to primitive type: int
    * org.mockito.exceptions.base.MockitoException:
    * Mockito cannot mock this class: class cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.
    * }}}
    *
    * An alternative workaround was to use `Mockito.doReturn(retVal).when(mockObj).someMethod`.
    */
  def mockWithDefaults[A: ClassTag]: A =
    Mockito.mock(
      classTag[A].runtimeClass.asInstanceOf[Class[A]],
      Mockito.withSettings().defaultAnswer(Mockito.RETURNS_DEFAULTS)
    )

  def capture[A: ClassTag]: ArgumentCaptor[A] =
    ArgumentCaptor.forClass(classTag[A].runtimeClass.asInstanceOf[Class[A]])
}

object MockSugar extends MockSugar
