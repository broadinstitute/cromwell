version draft-3

import "badly_named_struct.wdl" alias StructCollision as StructCollision2

struct A {
  Int i
  Float f
}

struct B {
  A a
}

struct StructCollision {
  Int i
}
