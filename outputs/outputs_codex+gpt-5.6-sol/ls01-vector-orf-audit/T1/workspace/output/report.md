# Construct ORF audit

The controlled T1 skills were freshly installed and opened in the specified order. They informed local DNA/ORF and synthetic-fusion checks; their optional remote tools were not called, and no absent vector features were inferred.

- **c01:** START, STOP, FRAME pass; TAG fails because terminal `TAA` stops translation before the downstream C-terminal FLAG tag.
- **c02:** START and TAG pass; STOP fails because the insert ends in `GAA`; FRAME fails because insert length is not divisible by three.
- **c03:** START and STOP pass; FRAME fails because the supplied claim is `out_of_frame`; TAG fails because terminal `TGA` prevents the C-terminal fusion.

All constructs fail overall, with issue labels respectively `TAG`, `STOP;FRAME`, and `FRAME;TAG`.
