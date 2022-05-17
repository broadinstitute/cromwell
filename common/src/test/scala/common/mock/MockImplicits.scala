package common.mock

import org.mockito.Mockito
import org.mockito.stubbing.OngoingStubbing

trait MockImplicits {

  /**
    * Ported from specs2's abandoned mock classes, provides DSL methods.
    *
    * https://github.com/etorreborre/specs2/commit/6d56660e70980b5958e6c4ed8fd4158bf1cecf70#diff-a2627f56c432e4bc37f36bc56e13852225813aa604918471b61ec2080462d722
    */
  implicit class MockEnhanced[A](methodCall: A) {
    def returns(result: A): OngoingStubbing[A] = {
      Mockito.when(methodCall).thenReturn(result)
    }

    def answers(function: Any => A): OngoingStubbing[A] = {
      Mockito.when(methodCall) thenAnswer {
        invocationOnMock => {
          val args = invocationOnMock.getArguments
          // The DSL behavior of the below is directly taken with thanks from the link above.
          args.size match {
            case 0 =>
              function match {
                case function0: Function0[_] =>
                  function0.apply().asInstanceOf[A]
                case _ =>
                  function.apply(invocationOnMock.getMock)
              }
            case 1 =>
              function(args(0))
            case _ =>
              function(args)
          }
        }
      }
    }

    def responds(f: Any => A): OngoingStubbing[A] = answers(f)
  }

}
