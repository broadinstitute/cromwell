# The 'volatile' Optimization

Available in Cromwell version 49 and higher.

### Effect on Call Caching:

The 'volatile' optimization is applied to tasks in their `meta` section.
Call caching will be disabled for any call to that task during the execution of the workflow. 

This is particularly useful if:

* One task can produce stochastic results but you still want to use call caching in the rest of the workflow.
* You want to guarantee that a task is never call cached for any other reason.

## Language Support

### WDL

In a WDL `task`, this optimization is specified by adding a `volatile` field to 
the task's `meta` section. Here's an example:

```wdl
version 1.0
 
task make_random_int {
  
  meta {
    volatile: true
  }
  
  command <<<
    echo $RANDOM
  >>>

  output {
    Int random = read_string(stdout())
  }
}
```

## Backend Support

The volatile keyword applies equally to all backends.
