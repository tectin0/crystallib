# bvparm2020.cif: not covered by crystallib's license

`bvparm2020.cif` is I. David Brown's accumulated table of bond valence parameters (version of
2020-11-25), distributed by the IUCr at
<https://www.iucr.org/resources/data/datasets/bond-valence-parameters>. It is included unchanged
(SHA-256 `ec073207e070f8e42e49734d11ebe93de4012cc0af07016a99920877688e9a68`) and is not part of
crystallib's GPL-3.0 license. Its own terms, from the file's header:

> This table may be used and distributed without fee for non-profit purposes providing
> 1) that this copyright notice is included and
> 2) no fee is charged for the table and
> 3) details of any changes made in this list by anyone other than the copyright owner are
> suitably noted in the _audit_update record.
> Please consult the copyright owner regarding any other uses.

crystallib compiles the table in only with the `bvparm2020` feature, which is off by default.
Turn it on only where these terms allow it, i.e. for non-profit use. Without the feature, crystallib
uses the bond valence parameters of O'Keeffe & Brese (1991) in `bond-valence-1991.toml`, which
cover every pair of elements less precisely.
