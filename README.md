# everything

## Book Blurb

This repository is meant to contain "Everything...or at Least Some of
It".
This is where I keep my notes from past coursework and literature
reviews.
These notes are written so that
[Tree of Knowledge (`tok`)](https://github.com/vgandari/tok) generates a
structured, textbook-like document from my notes.

## How to Generate the Book

The document is generated from a subset of the files in the `main`
directory containing text for some topics of interest.

Here's what happens when you run `tok`:

- `tok` builds a graph of all topics relevant to the topics of interest
  based on the dependency relationships declared in the files.
- `tok` sorts the topics in the "best" order while respecting dependency
 relationships declared in the files
- `tok` calls LaTeX to generate a textbook-like PDF document that can be
  used for reference.

The resulting document will include all, and only all relevant topics
leading up to the topic(s) of interest.
To include _everything_, run

```sh
tok $(find . -name '*.yaml' -print) -wu --headings
```

in the `main` directory.

For shorter documents, the `--extra-headings` option may make for nicer
chapter/section heading generation and placement.

See the documentation for
[Tree of Knowledge (tok)](https://github.com/vgandari/tok) for help and
more options.
