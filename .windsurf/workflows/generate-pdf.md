---
description: Generate and process a new PDF version of the TEP manuscript
---

Generate a new PDF version of the TEP-WB manuscript from the built site, with compression and metadata embedding.

## Quick Method (Universal Generator)
```bash
cd /Users/matthewsmawfield/www/TEP-WB
python /Users/matthewsmawfield/www/TEP/scripts/generate_pdf_universal.py . --quality maximum --wait-time 5
```

## Manual Steps

1. **Build the static site** - Generate the latest static HTML with all components:
   ```bash
   cd site && node build.js
   ```

2. **Generate PDF from site** - Convert the built HTML to PDF using Playwright (maximum quality):
   ```bash
   python scripts/utils/generate_site_pdf.py --quality maximum --wait-time 5
   ```

3. **Process PDF with metadata** - Compress and embed academic metadata (DOI, abstract, keywords):
   ```bash
   python scripts/utils/process_pdf.py site/public/docs/13-TEP-WB-v0.2-Kilifi.pdf --quality ebook
   ```

4. **Copy compressed PDF to root** - Copy the compressed PDF to project root:
   ```bash
   cp site/public/docs/13-TEP-WB-v0.2-Kilifi.pdf ./13-TEP-WB-v0.2-Kilifi.pdf
   ```

5. **Verify both PDFs** - Check both compressed PDFs exist and display their properties:
   ```bash
   ls -lh site/public/docs/13-TEP-WB-v0.2-Kilifi.pdf ./13-TEP-WB-v0.2-Kilifi.pdf
   exiftool -Title -Author -Creator site/public/docs/13-TEP-WB-v0.2-Kilifi.pdf 2>/dev/null | head -10
   ```
