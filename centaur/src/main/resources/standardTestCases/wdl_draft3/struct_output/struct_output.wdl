version 1.0

# This workflow has been arbitrarily selected to test that tabs work as whitespace.
# Be very careful while editing, as many editors automatically convert spaces to tabs.

struct FooStruct {
	Int simple
	Pair[Array[Int], Map[String, Boolean]] complex
}

workflow struct_output {
	output {
		FooStruct myFoo = object {
			simple: 5,
			complex: ([5], { "t": true })
		}
	}
}
