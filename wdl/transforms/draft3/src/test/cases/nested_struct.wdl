version 1.0

struct A {
  Int i
  Float f
}

struct B {
  A a
  Int i
  Float f
}

workflow nested_struct {
  B b = object {
    a: object { i: 5, f: 5.5 },
    i: 6,
    f: 6.6
  }

  output {
    Float f = b.a.f
  }
}
