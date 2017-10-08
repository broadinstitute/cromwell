# wdl4s wdl_refactor Branch

This was created based on 0.19 code with the goal to clean up the API and push WDL-specific code that existed in Cromwell down into the wdl4s library.

## Goals / Implementation Plan

### What's currently implemented

* Cleaner API and separate out lexical vs. graph relationships
* Nested scatter support
* If statement support
* Global variables support
* Much more comprehensive fully-qualified name definition with corresponding syntax errors
* Declarations as graph nodes -- calls can depend on calls indirectly through declarations
* Better support for operations like generating lookup functions (e.g. `Call.lookupFunction`) with scope resolution.
* More comprehensive test suite (1897 tests!)
* Refactored the test suite, created file-based test cases for graph and lexical relationships.
* 70% of change is debt reduction, 30% new features.

### Immediate plan:

* Things like lookup function generation and variable resolution need more tests and more design.  probably 70% there right now.  My goal was to push as much code related to expression lookup functions and WDL functions into wdl4s as possible.  Cromwell shouldn't need to do wdl-specific lexical rules to resolve a variable.  There are some subtleties around resolving variables in the context of a `Task` versus a `Call` that might need some more in depth tests.
* Figure out specific if statement semantics and update the WDL specification.  See `src/test/cases/if_statement/test.wdl` and `src/test/cases/scatter_within_if/test.wdl` for some starting points about the discussion of if-statement semantics.  Some things that need to be figured out:
** How do you handle a call that depends on another call that's within an if-statement?  Use optional inputs for this?
** What if there's an if-statement within a scatter block?  What happens to those shards that are NOT run?
* Discuss and decide on specifics of nested scatter semantics, document in WDL specification
* Integrate the changes above with Cromwell

### Longer-term plan

* Design / Implement sub-workflows with `call workflow_name` syntax.

## Lexical vs. Graph Relationships

The API is now much better at distinguishing graph and lexical relationships.

Graph relationships describe dependencies between nodes.

Lexical relationships describe relationships between elements in the source code.  For example, if a `call` statement is within a `scatter` block, then the call's *parent* is the scatter block.

```
workflow w {
  call x
  call y {input: i=x.o}
}
```

* `x` and `y` are **siblings** whose **parent** is `w`
* `y` is **downstream** of `x` and `x` is **upstream** of `y`

### Graph relationships

I introduced the concept of a `GraphNode` trait:

```scala
trait GraphNode {
  def upstream: Set[GraphNode]
  def downstream: Set[GraphNode]
}
```

The following classes are `GraphNode`s:

* `Call`
* `Declaration`
* `If`
* `Scatter`

`GraphNode` describes relationships to other `GraphNode`s.  `GraphNode.upstream` returns a set of `GraphNode`s that the current `GraphNode` depend on.  `GraphNode.downstream` returns the set of `GraphNode`s that depends on the current node.

For example

```
task a {
  command { ... }
  output { Int o = read_int(stdout()) }
}

task b {
  Array[Int] ints
  command { ... }
  output { Int o = read_int(stdout()) }
}

workflow w {
  call a as a0
  call a as a1
  Array[Int] ints = [a0.o, a1.o]
  call b {input: ints=ints}
}
```

`b.upstream` evaluates to a set of one element: the declaration for `ints` in the workflow.  `upstream` and `downstream` only return `GraphNode`s that they *directly depend on*.

If you'd like `b.upstream` to return `call a as a0`, `call a as a1`, AND `Array[Int] ints = [a0.o, a1.o]`, this can be done easily like this:

```scala
def upstreamRecursive(node: GraphNode): Set[GraphNode] = {
  node.upstream ++ node.upstream.flatMap(upstreamRecursive)
}
```

### Lexical Relationships

A lexical relationship describes how `Scope`s relate to each other.  This relationship is always tied to how they are represented in source code.  Currently the following classes are `Scope`s:

* `Call`
* `Declaration`
* `Task`
* `Workflow`
* `Namespace`
* `If`
* `Scatter`

Lexical operations:

```scala
trait Scope {
  def parent: Option[Scope]
  def children: Seq[Scope] // in source code order
  def ancestry: Seq[Scope]
  def descendants: Set[Scope]
  def closestCommonAncestor(other: Scope): Option[Scope]
  def resolveVariable(name: String): Option[GraphNode]
  def lookupFunction(inputs: WorkflowCoercedInputs,
                     wdlFunctions: WdlFunctions[WdlValue],
                     shards: Map[Scatter, Int] = Map.empty[Scatter, Int],
                     relativeTo: Scope = this): String => WdlValue
}
```



If statements
=============

```
task A {
  String i
  command {ps}
  output {String o = read_string(stdout())}
}

workflow w {
  Int i
  Array[String] arr

  call A

  if (i == 2) {
    call A as B
  }

  if (A.o == "foo") {
    call A as C
  }

  if (A.o == "bar") {
    scatter(x in arr) {
      call A as D {input: i=x}
    }
  }

  call A as E {input: i=C.o}
}
```

For the above WDL:

```scala
val d = namespace.calls.find(_.unqualifiedName == "D").get

println(d.ancestry.map(_.fullyQualifiedNameWithIndexScopes))
// List(w.$if_2.$scatter_0, w.$if_2, w, )

println(d.upstream.map(_.fullyQualifiedNameWithIndexScopes))
// Set(w.$if_2.$scatter_0)

println(upstreamRecursive(d).map(_.fullyQualifiedNameWithIndexScopes))
//Set(w.$if_2.$scatter_0, w.arr, w.$if_2, w.A)
```

Nested Scatters
===============

```
task inner {
  Int i
  command { echo ${i} }
  output {
    Int out = read_int(stdout()) + 1
  }
}

task outer {
  Array[Array[Int]] matrix
  command { cat ${write_tsv(matrix)} }
  output {
    String tsv = read_string(stdout())
  }
}

workflow w {
  Array[Array[Int]] array = [[0,1,2],[3,4,5],[6,7,8]]

  scatter(i in array) {
    scatter(j in i) {
      call inner {input: i=j}
    }
  }

  # inner.out should be [[1,2,3],[4,5,6],[7,8,9]] (???)
  call outer {input: matrix=inner.out}
}
```

Now doing some graph and lexical operations on the `call inner`:

```
val inner = namespace.calls.find(_.unqualifiedName == "inner").get

println(inner.ancestry.map(_.fullyQualifiedNameWithIndexScopes))
// List(w.$scatter_0.$scatter_1, w.$scatter_0, w, )

println(inner.upstream.map(_.fullyQualifiedNameWithIndexScopes))
// Set(w.$scatter_0.$scatter_1)

println(inner.downstream.map(_.fullyQualifiedNameWithIndexScopes))
// Set(w.outer)

println(inner.resolveVariable("j").map(_.fullyQualifiedNameWithIndexScopes))
// Some(w.$scatter_0.$scatter_1)

println(upstreamRecursive(inner).map(_.fullyQualifiedNameWithIndexScopes))
// Set(w.$scatter_0.$scatter_1, w.$scatter_0, w.array)
```

## Fully-Qualified Names

adder.wdl

```
task adder {
  Int i
  command { echo ${i} }
  output { Int o = read_int(stdout()) + 1 }
}

workflow add_workflow {
  call adder
}
```

main.wdl

```
import "adder.wdl" as ns

workflow main {
  call ns.adder
}
```

List of FQNs in `main.wdl`

|Fully-Qualified Name|
|--------------------|
|`main`|
|`main.adder`|
|`main.adder.i`|
|`ns`|
|`ns.adder`|
|`ns.adder.i`|
|`ns.add_workflow`|
|`ns.add_workflow.adder`|
|`ns.add_workflow.adder.i`|
|``|

> **TODO**: `main.adder.o` should also be a valid FQN

Also detects name collisions:

```
import "foobar" as ns
Int main = 2
workflow main {
  call ns.adder
}
```

Results in:

```
ERROR: Sibling nodes have conflicting names:

Declaration defined here (line 2, col 1):

Int main = 2
^

Workflow statement defined here (line 3, col 10):

workflow main {
         ^
```

## Global Variables

```
Int global = 2

task t {
  Int i = global + 2
  command { echo ${i} ${global} }
  runtime {
    memory: global + "GB"
  }
  output { String o = read_string(stdout()) }
}

workflow w {
  Int x = global + 2
  call t
}
```

|Fully-Qualified Name|
|--------------------|
|`global`|
|`t`|
|`t.i`|
|`w`|
|`w.x`|
|`w.t`|
|`w.t.i`|

## Test Suite Redesign

Tests for lexical relationships, graph relationships, and FQNs have been refactored to be auto-generating.

* `src/test/cases` contains 1 directory per WDL to test.  Each directory **must** contain a `test.wdl` file.
* `src/test/cases/*/*.json` are the expected results.  For example, `src/test/cases/nested_scatter/downstream.json` contains a JSON object where keys are FQNs in `./test.wdl` and values are array of downstream FQNs.

```
$ ls -lah src/test/cases
total 12K
drwxr-xr-x 14 sfrazer CHARLES\Domain Users 476 Aug  5 10:02 .
drwxr-xr-x  8 sfrazer CHARLES\Domain Users 272 Aug  2 08:54 ..
-rw-r--r--  1 sfrazer CHARLES\Domain Users 154 Aug  2 08:54 cgrep.wdl
drwxr-xr-x 11 sfrazer CHARLES\Domain Users 374 Aug  4 08:26 command_parameters
drwxr-xr-x 11 sfrazer CHARLES\Domain Users 374 Aug  4 12:46 depend_on_declaration
drwxr-xr-x 11 sfrazer CHARLES\Domain Users 374 Aug  2 08:54 if_statement
drwxr-xr-x 11 sfrazer CHARLES\Domain Users 374 Aug  2 08:54 nested_scatter
drwxr-xr-x 11 sfrazer CHARLES\Domain Users 374 Aug  5 12:39 nested_scatter2
-rw-r--r--  1 sfrazer CHARLES\Domain Users  75 Aug  2 08:54 ps.wdl
drwxr-xr-x 11 sfrazer CHARLES\Domain Users 374 Aug  4 10:24 scatter_within_if
drwxr-xr-x 11 sfrazer CHARLES\Domain Users 374 Aug  4 08:26 three_step
drwxr-xr-x 11 sfrazer CHARLES\Domain Users 374 Aug  4 08:26 three_step_import
drwxr-xr-x 11 sfrazer CHARLES\Domain Users 374 Aug  5 13:16 three_step_import_namespace
-rw-r--r--  1 sfrazer CHARLES\Domain Users 119 Aug  2 08:54 wc.wdl
```

```
$ ls -lah src/test/cases/nested_scatter
total 36K
drwxr-xr-x 11 sfrazer CHARLES\Domain Users  374 Aug  2 08:54 .
drwxr-xr-x 14 sfrazer CHARLES\Domain Users  476 Aug  5 10:02 ..
-rw-r--r--  1 sfrazer CHARLES\Domain Users  713 Aug  2 08:54 ancestry.json
-rw-r--r--  1 sfrazer CHARLES\Domain Users  526 Aug  2 08:54 children.json
-rw-r--r--  1 sfrazer CHARLES\Domain Users  436 Aug  2 08:54 downstream.json
-rw-r--r--  1 sfrazer CHARLES\Domain Users  516 Aug  2 08:54 fqn.json
-rw-r--r--  1 sfrazer CHARLES\Domain Users  648 Aug  2 08:54 fqn_index_scopes.json
-rw-r--r--  1 sfrazer CHARLES\Domain Users    4 Aug  2 08:54 inputs.json
-rw-r--r--  1 sfrazer CHARLES\Domain Users  474 Aug  2 08:54 parents.json
-rw-r--r--  1 sfrazer CHARLES\Domain Users 1.1K Aug  2 08:54 test.wdl
-rw-r--r--  1 sfrazer CHARLES\Domain Users  444 Aug  2 08:54 upstream.json
```

```
$ cat src/test/cases/nested_scatter/downstream.json
{
  "B.B_in": [],
  "C.C_in": [],
  "D.D_in": [],
  "w.$scatter_0": ["w.$scatter_1", "w.$scatter_2", "w.B", "w.C", "w.E"],
  "w.$scatter_1": ["w.G"],
  "w.$scatter_2": ["w.H"],
  "w.$scatter_3": ["w.F"],
  "w.A": ["w.$scatter_0", "w.$scatter_3"],
  "w.B": ["w.$scatter_1", "w.$scatter_2", "w.C", "w.D"],
  "w.B.B_in": [],
  "w.C": [],
  "w.C.C_in": [],
  "w.D": [],
  "w.D.D_in": [],
  "w.E": [],
  "w.F": [],
  "w.G": [],
  "w.H": []
}
```

How to use this:

* The goal is to have all WDL files that are used in tests to also automatically have their lexical relationships, graph relationships, inputs, and FQNs tested.
* If you want to add a WDL file, create a new directory in `src/test/cases`.  Then, in that directory create `test.wdl`.  Run `sbt test`.  The directory will now be populated with a bunch of `.json` files.  Verify that the contents of these is correct, then add them to Git.
