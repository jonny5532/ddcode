ddcode
======

2D Barcode Reader Library

This library provides routines for decoding PDF417 2D barcodes.

Currently implemented is a Reed-Solomon implementation suitable for generating and correcting PDF417 codes.

The Reed-Solomon decoder uses the Berlekamp-Massey/Fourney/Chien method, providing error correction for up to t errors (where there are 2t added parity codes). Erasure correction is currently not implemented.