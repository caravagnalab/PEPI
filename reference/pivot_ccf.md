# Stan data construction for the unified multirates model Pivot a long-format CCF table into an array shaped for Stan.

Stan data construction for the unified multirates model Pivot a
long-format CCF table into an array shaped for Stan.

## Usage

``` r
pivot_ccf(long_df, id_col, id_levels, time_levels, epistate = TRUE)
```

## Arguments

- long_df:

  Long-format tibble with an id column, \`time\`, \`ccf\` (and
  \`epistate\` if \`epistate = TRUE\`).

- id_col:

  Name of the id column in \`long_df\`.

- id_levels:

  Character vector of ids, in the order they should appear along the
  first array dimension.

- time_levels:

  Numeric vector of sampling times, in ascending chronological order
  (first = earliest sample).

- epistate:

  If TRUE, returns a \`\[length(id_levels), length(time_levels), 2\]\`
  array (epistate \`-\`/\`+\` on the third dimension); otherwise a
  \`\[length(id_levels), length(time_levels)\]\` matrix.

## Value

An array/matrix of CCF values.
