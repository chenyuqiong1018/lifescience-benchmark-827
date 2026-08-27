# Construct ORF audit

All three constructs were checked only against the supplied fields and frozen rule. No vector or sequence features were inferred.

- **c01:** start, stop, and frame checks pass. The terminal `TAA` would stop translation before the downstream C-terminal FLAG tag, so TAG fails. Overall: fail.
- **c02:** start and tag compatibility checks pass. The insert ends in `GAA`, not a stop codon, so STOP fails; its length is not divisible by three, so FRAME also fails despite the supplied `in_frame` claim. Overall: fail.
- **c03:** start and terminal-stop checks pass. The supplied claim is `out_of_frame`, so FRAME fails; the terminal `TGA` also prevents translation into the C-terminal FLAG tag, so TAG fails. Overall: fail.

Issue labels are emitted in the fixed order START, STOP, FRAME, TAG.
