bucket = 'gs://bt-343-testing/'
prefix = bucket + \
    'lorem-ipsum-dolor-sit-amet/consectetur-adipiscing-elit/nullam-in-aliquet-sapien/phasellus-at-feugiat-diam'


def main():
    output = open('/Users/mcovarr/inputs.tsv', 'w')

    # Write an array of arrays of input files.
    # The outer array should be ~700 wide, the inner array should be ~6 * 285 or about 1700.
    # Cycle through all of the inputs so a hash is requested for all of them (avoid the
    # root workflow file hash cache actor coalescing these requests).

    # scatter_width = 700
    # inputs_per_call = 1700
    scatter_width = 3
    inputs_per_call = 3

    lines = []
    for a in range(20):
        for b in range(10):
            for c in range(10):
                lines.append(f'{prefix}/{a}-{"a"*64}/{b}-{"b"*64}/{c}-{"c"*64}/input.txt')

    num_inputs = len(lines)

    for i in range(scatter_width):
        a = (inputs_per_call * i) % num_inputs
        b = (inputs_per_call * (i + 1)) % num_inputs

        if a < b:
            raw = lines[a:b]
        else:
            raw = lines[:b] + lines[a:]

        output.write('\t'.join(raw))
        if i < scatter_width - 1:
            output.write('\n')

    output.close()


if __name__ == '__main__':
    main()
