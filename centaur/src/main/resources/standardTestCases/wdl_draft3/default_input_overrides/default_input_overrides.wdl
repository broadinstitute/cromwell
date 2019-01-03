version 1.0

workflow default_input_overrides {

  meta {
    description: "Checks that workflow input overriding works as expected by https://github.com/openwdl/wdl/pull/161"
  }

  parameter_meta {
    a: "Is provided as an input in the inputs JSON (10101)"
    b: "Is provided as an input in the inputs JSON (20202)"
    c: "Not in the inputs JSON so picks up a default from an intermediate value (75)"
    d: "Not in the inputs JSON so picks up a default from an intermediate value (75)"

    a2: "Is provided as an input in the inputs JSON (30303)"
    b2: "Gets the provided input value of b as its default (20202)"
    c2: "Is provided as an input in the inputs JSON (40404)"
    d2: "Gets the default value of d as its default (75)"

    z: "No default so we rely on having a value from inputs (50505)"
    z2: "Gets the provided input value of z as its default via an intermediate value (50505)"
  }

  input {
    Int a = some_intermediate_value
    Int b = some_intermediate_value
    Int c = some_intermediate_value
    Int d = some_intermediate_value
    Int a2 = a
    Int b2 = b
    Int c2 = c
    Int d2 = d

    Int z
    Int z2 = intermediate_z
  }

  Int some_intermediate_value = 75
  Int intermediate_z = z

  output {
    Int a_out = a
    Int b_out = b
    Int c_out = c
    Int d_out = d

    Int a2_out = a2
    Int b2_out = b2
    Int c2_out = c2
    Int d2_out = d2

    Int z_out = z
    Int z2_out = z2
  }
}
