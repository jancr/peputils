# Peputils - Peptide Utility Library, 

Primarily a collection of Panda Extensions for Peptide Data analysis.

## How to use

By importing the library

```
import peputils
```

A series of `peptidomics` extensions methods are available for pandas peptidomics data

```
import pandas as pd
df = pd.DataFrame.peptidomics.load_upf_meta(upf_file, meta_file, "Mouse Brain")
annotations = df.peptidomics.create_annotation()
```

