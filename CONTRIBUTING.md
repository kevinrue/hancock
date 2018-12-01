
## Considerate contributions are more than welcome!

Here are a few guidelines to help develop and maintain a consistent coding style.

## Primary objectives

Most importantly, contributions should directly contribute to the primary objectives of the package. namely:

1. Apply signatures to assign cell identities to new data sets
2. Learn new signatures from data sets, in a format compatible with (1.)

## Internal functions

Internal functions should also be documented using roxygen comments (http://r-pkgs.had.co.nz/man.html).
However, those do not have to be as comprehensive as exported functions.
Nevertheless, required sections are:

- A title
- `@rdname INTERNAL_<...>` with `<...>` being the name of the function (without any trailing ".").
    Make sure that your `.gitignore` contains the entry `INTERNAL_*`. Do _not_ push INTERNAL documentation online.
- `@param`
- `@return`
- `@author`
