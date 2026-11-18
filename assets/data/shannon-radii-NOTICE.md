# shannon-radii.csv: from pymatgen, MIT license

`shannon-radii.csv` holds the effective ionic and crystal radii of R. D. Shannon, "Revised
effective ionic radii and systematic studies of interatomic distances in halides and
chalcogenides", Acta Crystallographica A32 (1976) 751–767, by element, oxidation state,
coordination (in Roman numerals; `SQ` square planar, `PY` pyramidal) and spin state.

It is pymatgen's `dev_scripts/periodic_table_resources/Shannon_Radii.csv`
(<https://github.com/materialsproject/pymatgen>, downloaded 2026-11-18), with its byte order mark
removed, its line endings changed to LF, and the charge 3 filled in on the row of Eu with
coordination VII, which lacked it (Eu³⁺, as in Shannon's table, between Eu³⁺ VI and VIII; Eu²⁺ VII
has its own row). The other rows are unchanged. pymatgen's license:

> The MIT License (MIT)
> Copyright (c) 2011-2012 MIT & The Regents of the University of California,
> through Lawrence Berkeley National Laboratory
>
> Permission is hereby granted, free of charge, to any person obtaining a copy of
> this software and associated documentation files (the "Software"), to deal in
> the Software without restriction, including without limitation the rights to
> use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
> the Software, and to permit persons to whom the Software is furnished to do so,
> subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
> FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
> COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
> IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
> CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
