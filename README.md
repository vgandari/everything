# everything

This repository is meant to contain everything...or at least some of it.
This is where I keep all the notes for everything I (am supposed to)
know.
Whenever I want to look up a topic, I check to see if I've written a
file for it.
If I have, I generate a PDF document that contains everything I've had
to learn (and may need to relearn) leading up to that topic or those
topics of interest.
I use my tool [Tree of Knowledge (tok)](https://github.com/vgandari/tok) to
generate a textbook-like document.
All topics in the document are guaranteed to be sotred in a way that
will make sense to the reader.

If you really want to include _everything_, run

```sh
tok $(find . -name '*.yaml' -print) -wu --headings
```

For shorter documents, the `--extra-headings` option may make for nicer
chapter/section heading generation and placement.

See the documentation for
[Tree of Knowledge (tok)](https://github.com/vgandari/tok) for help.
