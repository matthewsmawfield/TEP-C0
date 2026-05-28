---
description: Generate and process a new PDF version of the TEP-C0 manuscript
---

Generate a new PDF version of the TEP-C0 manuscript from the built site, with compression and metadata embedding.

## Quick Method (Site Generator)
```bash
python scripts/generate_site_pdf.py --quality maximum --wait-time 5
```

## Manual Steps

1. **Build the static site** - Generate the latest static HTML with all components:
   ```bash
   cd site && node build.js
   ```

2. **Generate PDF from site** - Convert the built HTML to PDF using Playwright (maximum quality):
   ```bash
   python scripts/generate_site_pdf.py --quality maximum --wait-time 5
   ```

3. **Process PDF with metadata** - Compress and embed academic metadata (DOI, abstract, keywords):
   ```bash
   python scripts/utils/process_pdf.py site/public/docs/14-TEP-C0-v0.1-Athens.pdf --quality ebook
   ```

4. **Copy compressed PDF to root** - Copy the compressed PDF to project root:
   ```bash
   cp site/public/docs/14-TEP-C0-v0.1-Athens.pdf ./14-TEP-C0-v0.1-Athens.pdf
   ```

5. **Verify both PDFs** - Check both compressed PDFs exist and display their properties:
   ```bash
   ls -lh site/public/docs/14-TEP-C0-v0.1-Athens.pdf ./14-TEP-C0-v0.1-Athens.pdf
   exiftool -Title -Author -Creator site/public/docs/14-TEP-C0-v0.1-Athens.pdf 2>/dev/null | head -10
   ```
