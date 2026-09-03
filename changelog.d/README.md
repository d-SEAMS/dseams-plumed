# News fragments

One file per user-visible change, in the same commit as the change.

Name: `+slug.<type>.md` or `<github-issue>.<type>.md`.

Types (Keep a Changelog): `added`, `changed`, `deprecated`, `removed`,
`fixed`, `security`.

```
towncrier create -c "DSEAMS_CAGES reports noccupied." +guests.added.md
```

The release cut is `towncrier build --version X.Y.Z`. That writes
`CHANGELOG.md` and removes the consumed fragments.
