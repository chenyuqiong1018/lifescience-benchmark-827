# Construct ORF audit

The selected whitelist skills informed local DNA validation, ORF/fusion compatibility checking, and reproducible computation. No external tools or sequence features were used.

- **c01:** START, STOP, and FRAME pass. Its terminal `TAA` prevents translation into the downstream C-terminal FLAG tag, so TAG fails. Overall: fail.
- **c02:** START and TAG pass. The insert ends in `GAA`, so STOP fails; its length is not divisible by three, so FRAME also fails. Overall: fail.
- **c03:** START and STOP pass. The supplied `out_of_frame` claim makes FRAME fail, and terminal `TGA` prevents translation into the C-terminal tag, making TAG fail. Overall: fail.

Every issue is emitted using only the frozen labels START, STOP, FRAME, and TAG.
