# Introduction

This specification aims to be a subset of the full [WDL specification](https://github.com/broadinstitute/common-workflow-language/blob/master/README.md) aimed at providing basic functionality to stakeholders.  Over time we plan on supporting more features from the WDL specification as they are needed or required.

For specifics on the grammar rules, refer to the original WDL specification.

# Overview

The WDL DSL is meant to be a *human readable and writable* way to express tasks and workflows.  The "Hello World" tool in CWL would look like this:

```
task hello {
  command {
    egrep '${pattern}' '${file in}'
  }
}
```

This describes a task, called 'hello', which has one parameter (message).  The value of the parameters (called "job parameters") is provided in a language specific way.  The reference implementation accepts the value of the parameters in JSON.  For example:

|Variable|Value    |
|--------|---------|
|pattern |^[a-z]+$ |
|in      |/file.txt|

Running the hello tool with these job parameters would yield a command line:

```
egrep '^[a-z]+$' '/file.txt'
```

A simple workflow that runs this task in parallel would look like this:

```
workflow example {
  array[file] files
  scatter(path in files) {
    call hello {input: in=path}
  }
}
```

## Task Definition

A task is purely declarative with a focus on constructing an atomic command which will be interpreted by the backend.  Tasks also specify their outputs and their corresponding types.  Tasks may also specify runtime and backend specific information.

### Sections

A task has exactly one `command` section, and optionally one `runtime` and/or `output` section

### Command Section

A command is a *task section* that starts with the keyword 'command', and encloses its body in curly braces or here-doc syntax.  The command specifies a template for a string that will be interpreted by the backend.  In a common case, this template would express a command line invocation, but is not required to.

### Command Parts

The parser should read characters from the command section (inside the curly braces or `<<<`) until it reaches a `${` character sequence.  This is interpreted as a literal string.

The parser should interpret anything enclosed in `${`...`}` as a variable.

Inside of the `${`...`}` must be declared a name of the variable (e.g. `${var}`).

A type may optionally be specified.  If one is not specified it defaults to `string` (e.g. `${file var}`, `${array[string] var}`)

The variable may have a quantifier suffix, which may be:

* `?` - means the variable is optional
* `*` - means the variable is optional or may accept a single element or an array of elements.
* `+` - means the variable may accept a single element or an array of elements

If the qualifier is `*` or `+`, then the `sep` attribute must be set.

Variables may specify attributes.  The only supported attribute currently is `sep` and must have a string value.

* `${sep=' ' file var+}` - may accept single file or `array[file]`.  If `array[file]` is given, combine elements with a space.
* `${'-f ' file var?}` - `var` is optional, if specified, prefix with `-f `.

```
grep '${start}...${end}' ${file input}
```

This command would be parsed as:

* `grep '` - string
* `${start}` - var
* `...` - string
* `${end}` - var
* `' ` - string
* `${file input}` - var

### Outputs Section

The outputs section defines which of the files and values should be exported after a successful run of this task.

Each output declares a variable with a corresponding type.

Each output must declare the value of that variable.

```
output {
  int threshold = read_int("threshold.txt")
  array[string] stuff = tsv("strings.tsv")
}
```

The task is expecting a URI called "threshold.txt".  Inside of that file must be one line that contains only an integer and whitespace.

This task also expects to read a URI "strings.tsv" and it expects it to contain a single-column TSV which will be interpreted as an array of strings (each line is an element in the array)

These functions, `read_int` and `tsv` and others, are implemented in an backend specific way.

### Runtime Section

The runtime section defines key/value pairs for runtime information needed for this task.  The key/value pairs are arbitrary and optional, from a language perspective.  However some backends may require certain runtime declarations.  Backends may interpret these runtime declarations in an implementation specific way.

For example, Herc may require memory, disk, and cores whereas another backend may not require any of these.

#### memory

Memory requirements for this task.  This should be a string value with an integer and suffixes like `B`, `KB`, `MB`, ... or binary suffixes `KiB`, `MiB`, ...

* "2MB"
* "6GiB"

#### disk

Disk requirements for the task.  Same format as `memory`.

#### cores

Floating point value for the number of cores or compute units that should be allocated to the task.

## Workflow Definition

A workflow is defined as the keyword `workflow` and the body being in curly braces.

```
workflow wf {
  array[file] files
  int threshold
  map[string, string] my_map
}
```

### Call Sub-Tasks

A workflow may call other tasks via the `call` keyword.  Following this keyword is the name of the task to call and optionally an alias if the name is already in use. 

```
workflow wf {
  call my_task
  call my_task as my_task_alias
  call my_task as my_task_alias2 {
    input: threshold=2
  }
}
```

The body of a `call` can declare variables and specify an `input` and `output` section.  The `input` section is key/value pairs where the key is the name of the variable inside the task and the value is the value to assign to that variable.  The `output` section works the same way, except for mapping outputs of the task to variables in the workflow.

As an example, here is a workflow in which the second task references an output from the first task:

```
task task1 {
  command {python do_stuff.py}
  output {file results = "stdout"}
}
task task2 {
  command {python do_stuff2.py ${file foobar}}
  output {file results = "stdout"}
}
workflow wf {
  call task1
  call task2 {input: foobar=task1.results}
}
```

### Scatter

A `scatter` clause has an iteration statement (e.g. `i in integers`) which describes which collection to scatter over and the name of each shard's element.  The body of the scatter clause declares which tasks and/or workflows to execute per shard.

```
scatter(i in integers) {
  call task1{input: num=i}
  call task2{input: num=task1.output}
}
```

In this example, `task2` depends on `task1`.  Tasks within a shard are completely self-contained therefore `task2` can only reference one possible `task1` for its inputs.

### Outputs

A workflow may contain an `output` section enclosed in curly braces.  Inside the body contains comma separated expressions for which outputs of tasks within the workflow to expose as outputs of the workflow.

```
workflow wf {
  call task1
  call task2
  output {
    task1.x, task2.y
  }
}
```

The outputs are defined as `task.output_variable`.

# Data Types

Primitives:
* string
* int
* float
* file
* uri
* boolean

Compound Types:
* array\[type\] (e.g. `array[int]`, `array[file]`)
* map\[type, type\] (e.g. `map[string, file]`)
