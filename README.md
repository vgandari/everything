# everything

## About

This is where I keep my notes.

My notes are distributed among a set of YAML files under the `yaml/`
directory.
Each YAML file contains LaTeX-formatted text and metadata to be
processed by [tok](https://github.com/vgandari/tok).

The resulting PDF is my book, _Everything...or at least some of it_.

## Requirements

This project requires XeLaTeX and
[tok](https://github.com/vgandari/tok).

## Building a PDF

Run the following for help on using `tok`.

```sh
tok --help
```

Once you run `tok` on a file or set of files, the output will be
`./yaml/output/main.pdf`.

If you would like to include all the notes in a comprehensive text,
you may do so by running, from the `yaml/` directory,

```sh
tok `find . -name '*.yaml' -print`
```

## Contributing

See the guidelines for how to best write content in a YAML file in the
documentation for [tok](https://github.com/vgandari/tok).

Pull requests are welcome from any domain of expertise.
Pull requests will be accepted if

- There is no plagiarism
- There are no errors when generating a textbook using
  [tok](https://github.com/vgandari/tok)
- There are no missing references in the generated pdf.

Please file an issue if, for example,

- Any of the above rules are broken
- Any information is inaccurate or unclear and you are not an expert in
  that domain (so that someone may submit a pull request with the
  correction)

Any feature requests for [tok](https://github.com/vgandari/tok) should
be submitted as an issue for that repository.
