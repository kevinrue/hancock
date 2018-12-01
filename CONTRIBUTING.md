
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
