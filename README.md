# emtrey (M-2PO, minimap 2 psl output)

[![Github release](https://img.shields.io/github/tag/rvolden/emtrey.svg?label=Version)](https://github.com/rvolden/emtrey/tags)

SAM to PSL conversion program. Preferably takes a SAM file produced by [minimap2](https://github.com/lh3/minimap2) to make a PSL file.

## Dependencies

- [Go](https://golang.org/dl/)

To compile, use `make`

### Usage

emtrey takes a minimap SAM file and outputs a PSL file.

Options:

```
-i    Input SAM file
-m    Use to indicate the use of a SAM file that came from minimap2
```

emtrey will output the PSL file to stdout.

Using a minimap sam file:

```bash
./emtrey -m -i minimap_output.sam >out.psl
```

Otherwise:

```bash
./emtrey -i alignment.sam >out.psl
```

--------------------------------------------------------------------------------

### Notes

- If you aren't using a minimap sam file, the first three columns of the psl file will probably be wrong, since most aligners use M to indicate matches and mismatches without using X or = in their CIGAR strings. However, matches and mismatches can be calculated using the `NM:i:` and `nn:i:` tags that minimap produces.

- The strand is determined using the 0x10 bit, as per the SAM file format spec. If you have a spliced minimap alignment, it will check for the extra tag `ts:A:`, which overwrites what was calculated from the bit.

- Be sure to leave the header in the SAM file, as it's needed for the target sizes.
