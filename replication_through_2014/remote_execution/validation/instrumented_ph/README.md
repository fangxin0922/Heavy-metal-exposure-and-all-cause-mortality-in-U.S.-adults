# Checked Peng–Huang native source

`src/ph_diagnostic.f` is the audited 33-field source used by the EMM analysis. It derives from the `crqfnb.f` and `rqfnb.f` sources in **quantreg 6.1**, with privately renamed routines and recorded numerical state. The R adapter checks the solver path and retains only the valid prefix.

Upstream sources:

- [quantreg 6.1 crqfnb.f](https://github.com/cran/quantreg/blob/6.1/src/crqfnb.f)
- [quantreg 6.1 rqfnb.f](https://github.com/cran/quantreg/blob/6.1/src/rqfnb.f)
- [quantreg 6.1 package metadata](https://github.com/cran/quantreg/blob/6.1/DESCRIPTION)

The package metadata identifies Roger Koenker as author/maintainer and Stephen Portnoy as a contributor to censored QR code, alongside the package's other contributors. Its license is **GNU GPL version 2 or later**. The retained `Output from Public domain Ratfor` comments identify the translator; they do not relabel the package-derived source as MIT. See the repository's GPL notice and license text.

Compile the supplied source directly on Linux, inside `src`:

```sh
R CMD SHLIB ph_diagnostic.f -o ph_diagnostic.so -llapack -lblas
```

The historical instrumentation generator is not included because it predates the final 33-field source. No compiled binary is distributed. Source MD5 from the completed run: `039fdc854e9e4489aff9277962d617a1`.

Source comments, native epsilon and iteration cap are preserved. Compilation and numeric-library versions can affect floating-point behavior. The full EMM workflow retains its documented solver-sensitivity qualifications.
