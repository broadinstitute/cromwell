# Little tool that uses a prime sieve to generate primes up to 100

cwlVersion: v1.0
$graph:
- id: iwdr_input_string_function
  class: CommandLineTool
  baseCommand: ['python3']
  requirements:
      - class: DockerRequirement
        dockerPull: "python:3.5.0"
      - class: InitialWorkDirRequirement
        listing:
            - entryname: ${var result = inputs.script_name + ".py"; return result;}
              entry: ${var scr = inputs.script; return scr + "print(result);"}

  stdout: "primes"
  inputs:
      - id: script
        type: string
        default: |
          import sys
          import math

          limit = int(sys.argv[1])
          sieve = [True for i in range(limit)]
          for i in range(2, math.floor(limit / 2)):
              if sieve[i]:
                  for j in range(i * 2, limit, i):
                      sieve[j] = False

          result = "["
          for i in range(2, limit):
              if sieve[i]:
                  if result != "[":
                      result += ", "
                  result += str(i)
          result += "]"

      - id: script_name
        type: string
        default: "prime_sieve.py"
  arguments: [ '${var result = inputs.script_name + ".py"; return result;}', '100' ]
  outputs:
      - id: prime_list
        type: string
        outputBinding:
          glob: primes
          loadContents: true
          outputEval: $(self[0].contents.trim())
