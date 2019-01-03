version 1.0

##
## A set of expressions which would cause errors, except luckily they're short-circuited
##
workflow draft3_short_circuit {

  if (false && 1/0 == 0) {
    Int a = 1
  }

  if (true || 1/0 == 0) {
    Int b = 1
  }

  Int c = if true then 1 else 1/0
  Int d = if false then 1/0 else 1

  output {
    Int? a_out = a
    Int? b_out = b
    Int c_out = c
    Int d_out = d
  }
}
