# Retry with More Memory

With this feature one can specify an array of strings which when encountered in the `stderr` file by Cromwell, 
allows the task to be retried with more memory. The retry will be counted against the `maxRetries` count mentioned in 
the `runtimeAtrributes` in the task. There are 2 settings for this feature:

* `system.memory-retry-error-keys` : the error keys that need to be set in Cromwell config
* `memory_retry_multiplier` : [optional] the factor by which the memory should be multiplied while retrying. This needs 
to be passed in through workflow options and should be in the range `1.0 ≤ multiplier ≤ 99.0` (note: if set to `1.0` the task
will retry with same amount of memory). If this is not specified, Cromwell will retry the task (if applicable) but not 
change the memory amount.

For example, if the error keys set in Cromwell config are as below and the multipler passed through workflow options is 
`"memory_retry_multiplier": 1.1` 
```hocon
system {
  memory-retry-error-keys = ["OutOfMemory", "Killed"]
}
```  
this tells Cromwell to retry the task with 1.1x memory when it sees either `OutOfMemoryError` or `Killed` in the `stderr` 
file. 

If the task has runtime attributes as below 
```hocon
runtimeAtrributes {
  memory: "1 GB"
  continueOnReturnCode: true
  maxRetries: 1
}
``` 
the task will be retried 1 more time if it runs out of memory, and this time with "1.1 GB". 

Please note that this feature currently only works in Google Cloud backend. Also, Pipelines API might adjust the 
memory value based on their standards for memory for a VM. So it's possible that even though the request says 1.1 GB 
memory, it actually allocated a bit more memory to the VM.

Two environment variables called `${MEM_UNIT}` and `${MEM_SIZE}` are also available inside the command block of a task,
making it easy to retrieve the new value of memory on the machine.
