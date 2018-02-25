# Little tool that uses a prime sieve to generate primes up to 100

cwlVersion: v1.0
$graph:
- id: iwdr_input_string
  class: CommandLineTool
  baseCommand: ['python3', 'prime_sieve.py', '100']
  requirements:
      - class: DockerRequirement
        dockerPull: "python:3.5.0"
      - class: InitialWorkDirRequirement
        listing:
            - entryname: $(inputs.script_name)
              entry: $(inputs.script)

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

          print(result)
      - id: script_name
        type: string
        default: "prime_sieve.py"
  arguments: [ $(inputs.script_name), '100' ]
  outputs:
      - id: prime_list
        type: string
        outputBinding:
          glob: primes
          loadContents: true
          outputEval: $(self[0].contents.trim())
