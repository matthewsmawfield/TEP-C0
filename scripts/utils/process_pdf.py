#!/usr/bin/env python3
"""Unified PDF Processing Script
Compresses PDF and embeds comprehensive metadata in one operation.

This script processes the TEP-C0 manuscript PDF ("Temporal Equivalence 
Principle: A Covariant Alternative to Cosmic Expansion") by compressing it 
for web distribution and embedding complete academic metadata for proper indexing and citation.

Usage:
    python process_pdf.py <input_pdf> [--quality ebook|printer|prepress|default]
    
Example:
    python process_pdf.py site/public/docs/14-TEP-C0-v0.1-Athens.pdf --quality ebook
"""

import subprocess
import sys
import os
from pathlib import Path
import argparse
import tempfile

from compress_pdf import compress_pdf as _compress_pdf


def compress_pdf(input_path, output_path, quality='ebook'):
    """Compress PDF using Ghostscript (wrapper with legacy return format)."""
    result = _compress_pdf(input_path, output_path, quality=quality)
    return {
        'original_mb': result['original_size'] / (1024 * 1024),
        'compressed_mb': result['compressed_size'] / (1024 * 1024),
        'reduction_pct': result['reduction_percent']
    }


def embed_metadata(pdf_path, metadata):
    """Embed metadata into PDF using exiftool."""
    cmd = ['exiftool']

    # Add all metadata fields
    for key, value in metadata.items():
        cmd.extend([f'-{key}={value}'])

    # Overwrite original
    cmd.extend(['-overwrite_original', pdf_path])

    try:
        subprocess.run(cmd, check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Exiftool metadata embedding failed: {e.stderr.decode()}")


def verify_metadata(pdf_path, expected_fields):
    """Verify metadata was embedded correctly."""
    cmd = ['exiftool'] + [f'-{field}' for field in expected_fields] + [pdf_path]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError:
        return None


def main():
    parser = argparse.ArgumentParser(
        description='Compress PDF and embed metadata in one operation'
    )
    parser.add_argument('input_pdf', help='Path to input PDF file')
    parser.add_argument(
        '--quality',
        choices=['screen', 'ebook', 'printer', 'prepress', 'default'],
        default='ebook',
        help='Compression quality (default: ebook)'
    )
    parser.add_argument(
        '--doi',
        default='10.5281/zenodo.20370144',
        help='DOI to embed in metadata'
    )

    args = parser.parse_args()

    input_path = Path(args.input_pdf).resolve()

    if not input_path.exists():
        print(f"Error: File not found: {input_path}")
        sys.exit(1)

    print(f"Processing: {input_path}")
    print(f"Quality: {args.quality}")
    print()

    # Step 1: Compress PDF
    print("Step 1: Compressing PDF...")
    with tempfile.NamedTemporaryFile(suffix='.pdf', delete=False) as tmp:
        tmp_path = tmp.name

    try:
        stats = compress_pdf(str(input_path), tmp_path, args.quality)

        # Replace original with compressed version
        os.replace(tmp_path, str(input_path))

        print(f"  Original:    {stats['original_mb']:.2f} MB")
        print(f"  Compressed:  {stats['compressed_mb']:.2f} MB")
        print(f"  Reduction:   {stats['reduction_pct']:.1f}%")
        print()

    except Exception as e:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        print(f"Error during compression: {e}")
        sys.exit(1)

    # Step 2: Embed metadata
    print("Step 2: Embedding metadata...")

    # Paper metadata - must match manuscript, CITATION.cff, and zenodo.txt
    metadata = {
        # Core identification
        'Title': 'Temporal Equivalence Principle: A Covariant Alternative to Cosmic Expansion',
        'Author': 'Matthew Lukin Smawfield',
        'Creator': 'Matthew Lukin Smawfield',

        # Scientific abstract with key results
        'Subject': (
            'This paper develops the Level-3 cosmological extension of the Temporal Equivalence '
            'Principle (TEP): the hypothesis that observational evidence normally interpreted as '
            'cosmic expansion may involve large-scale Temporal Shear. Utilizing a global Bayesian '
            'synthesis of 1,701 Pantheon+ supernovae and Planck 2018 acoustic anchors, the analysis '
            'finds that the Universal TEP model achieves decisive statistical superiority over '
            'standard LambdaCDM. High-fidelity nested sampling yields strong Bayesian evidence '
            '(Bayes Factor = 35.6) in favor of non-integrable temporal transport. By tracing the '
            'continuous screening transition across the cosmic web, the framework formally resolves '
            'the Hubble tension and independently validates structure growth against Planck sigma_8 '
            'measurements.'
        ),

        # Keywords for indexing
        'Keywords': (
            'Temporal Equivalence Principle; TEP; Temporal Cosmology; '
            'Cosmic Redshift; FLRW; Proper-Time Transport; Dark Energy; '
            'Supernovae; Bayesian Inference; Modified Gravity; Temporal Shear'
        ),

        # Production metadata
        'Producer': 'TEP-C0 Research Project - Version 0.1 (Athens)',

        # Rights and identifiers
        'Copyright': 'Creative Commons Attribution 4.0 International License (CC BY 4.0)',

        # Dates
        'CreationDate': '2026:05:25 00:00:00',
        'ModifyDate': '2026:05:25 00:00:00',

        # XMP Dublin Core metadata (exiftool uses these prefixes)
        'XMP-dc:Creator': 'Matthew Lukin Smawfield',
        'XMP-dc:Title': 'Temporal Equivalence Principle: A Covariant Alternative to Cosmic Expansion',
        'XMP-dc:Description': 'TEP cosmological extension: cosmic redshift as emergent proper-time transport',
        'XMP-dc:Rights': 'CC BY 4.0',
        'XMP-dc:Identifier': f'doi:{args.doi}',
        'XMP-dc:Source': 'https://github.com/matthewsmawfield/TEP-C0',
        'XMP-dc:Publisher': 'Zenodo',
        'XMP-dc:Date': '2026-05-25',
        'XMP-dc:Type': 'Preprint',
        'XMP-dc:Format': 'application/pdf',
        'XMP-dc:Language': 'en',

        # PRISM (Publishing Requirements for Industry Standard Metadata)
        'XMP-prism:DOI': args.doi,
        'XMP-prism:URL': 'https://github.com/matthewsmawfield/TEP-C0',
        'XMP-prism:VersionIdentifier': '0.1',
        'XMP-prism:PublicationName': 'TEP Research Series',

        # PDF/A metadata
        'XMP-pdfaid:Part': '1',
        'XMP-pdfaid:Conformance': 'B'
    }

    try:
        embed_metadata(str(input_path), metadata)
        print("  Metadata embedded successfully")
        print()

    except Exception as e:
        print(f"Error during metadata embedding: {e}")
        sys.exit(1)

    # Step 3: Verify
    print("Step 3: Verifying metadata...")
    verification = verify_metadata(
        str(input_path),
        ['Title', 'Author', 'Subject', 'Keywords', 'Creator', 'Copyright']
    )

    if verification:
        print("  ✓ Metadata verified")
        print()
        print("Verification output:")
        print(verification)
    else:
        print("  ⚠ Could not verify metadata")

    print()
    print(f"✓ Processing complete: {input_path}")
    print(f"  Final size: {os.path.getsize(input_path) / (1024 * 1024):.2f} MB")


if __name__ == '__main__':
    main()
