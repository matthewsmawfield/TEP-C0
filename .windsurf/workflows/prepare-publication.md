---
description: Prepare TEP-WB project for new publication release
---

# Prepare Publication Workflow

Comprehensive workflow to prepare the TEP-WB project for a new publication release. This ensures version consistency, date synchronization, abstract alignment, PDF generation, and final verification across the entire codebase.

## Critical Date Rules

- **First Published Date**: **NEVER CHANGES** - 19 March 2026 for TEP-WB (found in site/index.html)
- **Citation/Publication Dates**: Use first published date (2026-03-19)
- **Last Updated Date**: Always set to TODAY when publishing
- **Modified Time Meta**: Use today's date in ISO format (e.g., {TODAY_ISO})

## Phase 1: Extract Source of Truth Data

First, extract the authoritative version, codename, dates, and abstract from `site/index.html`:

```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB"

# Extract version and codename
echo "=== VERSION INFO ==="
grep -E "(Version:|version)" site/index.html | head -5

# Extract dates
echo "=== DATES ==="
grep -E "(First published|Last updated|published_time|modified_time)" site/index.html | head -10

# Extract first published date (for citations)
echo "=== FIRST PUBLISHED DATE (for citations) ==="
grep -E "article:published_time|citation_date|citation_publication_date" site/index.html | head -5
```

**Record these values:**
- Current Version: __________ (e.g., v0.5)
- Codename: __________ (e.g., Tortola)
- First Published: __________ (e.g., 19 December 2025 / 2026-03-19)
- Today's Date: __________ (e.g., {TODAY} / {TODAY_ISO})

## Phase 2: Update Version.txt with Git Changes

Update `version.txt` with a concise summary of changes since the last version:

```bash
# View git status to see what changed
git status

# View diff summary for key files
git diff --stat HEAD~1

# For version.txt, summarize in short bullet points:
# - Maximum 6 bullets (0-6 depending on changes)
# - High-level summary of latest changes since last version only
# - Focus on pipeline and manuscript changes
# - Exclude SEO, trivial formatting, or infrastructure changes
```

**Manually edit `version.txt`:**
```
v{X.Y} ({Codename}) - {Brief Summary}

- {Key pipeline/manuscript change}
- {Key pipeline/manuscript change}
# 0-6 bullets max, high-level summary of latest changes only
# Focus on pipeline and manuscript changes
# Exclude: SEO, trivial formatting, infrastructure
```

## Phase 3: Synchronize Version/Codename Across All Files

Update the following files with consistent version and codename:

### 3.1 VERSION.json (Root)
```json
{
  "version": "X.Y",
  "codename": "Codename",
  "description": "White Dwarf Cooling and the Temporal Equivalence Principle"
}
```

### 3.2 version.txt (Root)
Already updated in Phase 2. Verify it shows `vX.Y (Codename)`.

### 3.3 README.md
- Line 9: Version badge → `vX.Y (Codename)`
- Line 10: Date → `First published: 19 March 2026 · Last updated: DD MMMM YYYY`
- Line 151: BibTeX note → `Preprint vX.Y (Codename)`

### 3.4 site/index.html
Multiple locations to update:
- Line ~53: PDF alternate link → `/public/docs/13-TEP-WB-vX.Y-Codename.pdf`
- Line ~57: PDF href → `/public/docs/13-TEP-WB-vX.Y-Codename.pdf`
- Line ~86: `article:published_time` → Keep as `2026-03-19T00:00:00Z` (first published)
- Line ~89: `article:modified_time` → Update to today `{TODAY_ISO}T00:00:00Z`
- Line ~90: `og:updated_time` → Update to today
- Line ~138: `citation_date` → Keep as `2026-03-19` (first published)
- Line ~147: `citation_publication_date` → Keep as `2026/03/19` (first published)
- Line ~151: `citation_pdf_url` → `/public/docs/13-TEP-WB-vX.Y-Codename.pdf`
- Line ~947: Header version → `Version: vX.Y (Codename)`
- Line ~949: Dates → `First published: 19 March 2026 · Last updated: DD MMMM YYYY`

### 3.5 site/manifest.json
- Line 5: version → `vX.Y (Codename)`
- Line 7: last_updated → `DD MMMM YYYY`

### 3.6 site/package.json
- Line 3: version → `X.Y.0` (semver format)
- Line 5: comment → `Package version matches manuscript version vX.Y (Codename)`

### 3.7 site/CITATION.cff
- Line 12: version → `vX.Y (Codename)`
- Line 13: date-released → Keep as `2025-12-19` (first published)
- Line 29: preferred-citation title → `(Codename vX.Y)`

### 3.8 site/citation.json
- Line 4: title → Must match index.html exactly
- Line 17: note → `Preprint, Version X.Y (Codename)`

### 3.9 site/CITATION.bib
- Line 3: title → Must match index.html exactly
- Line 9: note → `Preprint, Version X.Y (Codename)`

### 3.10 site/codemeta.json
- Line 16: dateModified → Today's date `{TODAY_ISO}`
- Line 17: version → `X.Y`

## Phase 3.11: Check/Cleanup Old PDF Versions

```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB"

echo "=== Checking for old PDF versions ==="
ls -la *.pdf 2>/dev/null || echo "No PDFs in root"
ls -la site/public/docs/*.pdf 2>/dev/null || echo "No PDFs in site/public/docs"

# Note: PDF naming uses dash: 13-TEP-WB-v0.5-Tortola.pdf
# If version changed, remove old versions to avoid confusion
```

## Phase 3.5: Synchronize Title Across All Files

The paper title in `site/index.html` is the source of truth. Extract and verify consistency:

```bash
# Extract canonical title from index.html
echo "=== Canonical Title from site/index.html ==="
grep -E '<title>|og:title|citation_title|dc.title' site/index.html | head -5

# Check title in key files
echo ""
echo "=== Title in README.md ==="
grep -E '^# |^title:' README.md | head -2

echo ""
echo "=== Title in CITATION.cff ==="
grep '^title:' site/CITATION.cff

echo ""
echo "=== Title in codemeta.json ==="
grep '"name":' site/codemeta.json | head -1

echo ""
echo "=== Title in citation.json ==="
grep '"title":' site/citation.json

echo ""
echo "=== Title in manifest.json ==="
grep '"title":' site/manifest.json

echo ""
echo "=== Title in CITATION.bib ==="
grep 'title=' site/CITATION.bib
```

**Update title in these files to match site/index.html:**

**README.md** (Line 1):
- Title in main heading: `# Title Here`

**site/CITATION.cff** (Line 2):
- `title: "Exact Title from index.html"`
- Line 29: preferred-citation title → `(Codename vX.Y)` suffix

**site/citation.json** (Line 4):
- `"title": "Exact Title from index.html"`

**site/codemeta.json** (Line 5):
- `"name": "Exact Title from index.html"`

**site/manifest.json** (Line 2):
- `"title": "Exact Title from index.html"`

**site/CITATION.bib**:
- `title = {Exact Title from index.html}`

**site/index.html** (already source of truth, but verify):
- Line ~14-17: `<title>` tag
- Line ~22-23: `meta name="title"`
- Line ~64-66: `og:title`
- Line ~99-101: `twitter:title`
- Line ~134-136: `citation_title`
- Line ~183-185: `dc.title`

## Phase 4: Synchronize Abstract Across All Files

The source of truth is `site/components/abstract.html`. Extract and propagate to:

### 4.1 Extract Canonical Abstract
```bash
# View the abstract from components
cat site/components/abstract.html
```

Copy the text between `<p>` and `</p>` (excluding keywords line).

### 4.2 Update Abstract in These Files

**README.md** (Line 17):
- Replace content after `## Abstract` heading

**zenodo.txt** (Line 1):
- Must follow this exact template structure:
  ```
  {Abstract text from site/components/abstract.html}

  Website: https://mlsmawfield.com/tep/{project-shortname}
  Code Availability: https://github.com/matthewsmawfield/{REPO_NAME}

  Keywords: {relevant keywords matching the paper}

  Open Science Statement: This work is a preprint and is open to community review, ideas, and collaboration. All materials required for full reproducibility—including data downloads, analysis scripts, code, and manuscripts—are open-source. Feedback and contributions to further test these results are welcome.
  ```
- Example from TEP-JWST:
  - Website: https://mlsmawfield.com/tep/jwst
  - Code: https://github.com/matthewsmawfield/TEP-JWST
  - Keywords: Cosmology: early universe – Galaxies: high-redshift – etc.
- For TEP-WB: website is `https://mlsmawfield.com/tep/wb`, code is `https://github.com/matthewsmawfield/TEP-WB`
- Keywords should match those in CITATION.cff and index.html meta tags

**site/CITATION.cff** (Line 9-10):
- Update abstract field (YAML folded style with `>`)

**site/codemeta.json** (Line 6):
- Update description field (escaped JSON string)

**site/manifest.json** (Line 8):
- Update description field (shorter summary acceptable)

**site/index.html**:
- Line ~26: meta description
- Line ~68: og:description
- Line ~104: twitter:description
- Line ~192: dc:description

## Phase 5: Run PDF Generation Workflow

// turbo
```bash
# Build the static site first
cd site && node build.js

# Generate PDF with maximum quality
python scripts/generate_site_pdf.py --quality maximum --wait-time 5
```

## Phase 6: Verify PDF Output

### 6.1 Check PDF Locations
```bash
ls -lh "13-TEP-WB-vX.Y-Codename.pdf"
ls -lh site/public/docs/"13-TEP-WB-vX.Y-Codename.pdf"
```

### 6.2 Verify PDF Metadata
```bash
exiftool -Title -Author -Creator -Subject -Keywords "13-TEP-WB-vX.Y-Codename.pdf"
```

**Expected metadata:**
- Title: "White Dwarf Cooling and the Temporal Equivalence Principle"
- Author: "Matthew Lukin Smawfield"
- Keywords should include: physics, cosmology, gravitational lensing, dark matter, TEP-WB

### 6.3 Verify PDF Links
Ensure the PDF URL in these locations points to `site/public/docs/`:
- `site/index.html` line 53, 57, 151 (citation_pdf_url)
- Any download links in components

## Phase 7: Final Consistency Checks

### 7.1 Version String Scan
```bash
# Search for any remaining old version strings
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB"
echo "=== Checking for inconsistent version strings ==="
grep -r "v0\.[0-9]" --include="*.html" --include="*.md" --include="*.json" --include="*.txt" --include="*.bib" --include="*.cff" . 2>/dev/null | grep -v ".git" | grep -v node_modules | grep -v "vX.Y"

# Check for PDF filename patterns
grep -r "13-TEP-WB-v" --include="*.html" --include="*.md" . 2>/dev/null | grep -v ".git"

# Check citation file versions specifically
echo ""
echo "=== Version in citation files (should all match) ==="
echo "CITATION.bib:"
grep "note.*Version" site/CITATION.bib
echo "citation.json:"
grep '"note":' site/citation.json
echo "CITATION.cff:"
grep '^version:' site/CITATION.cff
echo "manifest.json:"
grep '"version":' site/manifest.json
```

### 7.2 Date Consistency Check
```bash
echo "=== First Published Dates (should all be Dec 2025) ==="
grep -r "2025-12-19\|19 December 2025" --include="*.html" --include="*.json" --include="*.md" --include="*.txt" . 2>/dev/null | grep -v ".git"

echo "=== Last Updated Dates (should be TODAY) ==="
grep -r "{TODAY}\|{TODAY_ISO}" --include="*.html" --include="*.json" --include="*.md" . 2>/dev/null | grep -v ".git"
```

### 7.3 Abstract Consistency Check
```bash
# Check abstract length in key files
echo "=== Abstract lengths ==="
echo "README.md:"
grep -A 1 "^## Abstract$" README.md | tail -1 | wc -c

echo "zenodo.txt:"
head -1 zenodo.txt | wc -c

echo "CITATION.cff:"
grep -A 5 "^abstract:" site/CITATION.cff | wc -c
```

### 7.4 zenodo.txt Template Check
```bash
echo "=== Checking zenodo.txt template structure ==="
echo "Website URL:"
grep -E "^Website:" zenodo.txt
echo "Code Availability:"
grep -E "^Code Availability:" zenodo.txt
echo "Keywords:"
grep -E "^Keywords:" zenodo.txt
echo "Open Science Statement:"
grep -E "^Open Science Statement:" zenodo.txt
echo ""
echo "Expected for TEP-WB:"
echo "  Website: https://mlsmawfield.com/tep/wb"
echo "  Code Availability: https://github.com/matthewsmawfield/TEP-WB"
echo "  Open Science Statement: Present"
```

### 7.5 Title Consistency Check
```bash
echo "=== Checking title consistency across files ==="

# Extract title from index.html (remove extra whitespace)
INDEX_TITLE=$(grep -oE '<title>.*</title>' site/index.html | sed 's/<title>//;s/<\/title>//' | tr -d '\n' | sed 's/  */ /g')
echo "Canonical title from index.html:"
echo "  $INDEX_TITLE"
echo ""

echo "Comparing with other files:"
echo ""
echo "README.md heading:"
grep -E '^# [^#]' README.md | head -1

echo ""
echo "CITATION.cff:"
grep '^title:' site/CITATION.cff

echo ""
echo "citation.json:"
grep '"title":' site/citation.json

echo ""
echo "codemeta.json:"
grep '"name":' site/codemeta.json | head -1

echo ""
echo "CITATION.bib:"
grep 'title=' site/CITATION.bib

echo ""
echo "WARNING: If any titles don't match exactly, update them to match index.html"
```

## Phase 8: Generate Commit Message

// turbo
```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB"

echo "========================================"
echo "SUGGESTED COMMIT MESSAGE"
echo "========================================"
echo ""
echo "vX.Y (Codename) - $(head -1 version.txt | sed 's/.*- //')"
echo ""
echo "Summary of Changes:"
tail -n +2 version.txt | grep "^- " || echo "- Version bump to vX.Y"
echo ""
echo "# Note: 0-6 bullets max, high-level summary of latest changes since last version"
echo "$(git diff --stat HEAD~1)"
echo ""
echo "========================================"
```

## Pre-Publication Checklist

### Version & Codename
- [ ] `VERSION.json` - version, codename correct
- [ ] `version.txt` - version header + changelog updated
- [ ] `README.md` - version, date updated
- [ ] `site/index.html` - version, PDF links, meta dates correct
- [ ] `site/manifest.json` - version, last_updated correct
- [ ] `site/package.json` - version, comment correct
- [ ] `site/CITATION.cff` - version correct, date-released is first published
- [ ] `site/citation.json` - title matches, note field version correct
- [ ] `site/CITATION.bib` - title matches, note field version correct
- [ ] `site/codemeta.json` - version, dateModified correct

### Dates
- [ ] First published date: 19 December 2025 (unchanged everywhere)
- [ ] Citation dates use first published date (2026-03-19/19)
- [ ] Last updated date: TODAY reflected in all locations
- [ ] Modified meta tags: TODAY in ISO format

### Title
- [ ] `site/index.html` - source of truth title correct
- [ ] `README.md` - title matches index.html
- [ ] `site/CITATION.cff` - title matches index.html
- [ ] `site/citation.json` - title matches index.html
- [ ] `site/codemeta.json` - name matches index.html
- [ ] `site/manifest.json` - title matches index.html
- [ ] `site/CITATION.bib` - title matches index.html

### Abstract
- [ ] `site/components/abstract.html` - source of truth current
- [ ] `README.md` - abstract matches
- [ ] `zenodo.txt` - abstract matches **AND** follows complete template (Website, Code, Keywords, Open Science Statement)
- [ ] `site/CITATION.cff` - abstract matches
- [ ] `site/codemeta.json` - description matches
- [ ] `site/manifest.json` - description matches
- [ ] `site/index.html` - meta descriptions match

### PDF
- [ ] PDF generated at root: `13-TEP-WB-vX.Y-Codename.pdf`
- [ ] PDF copied to site: `site/public/docs/13-TEP-WB-vX.Y-Codename.pdf`
- [ ] Old PDF versions removed from root and site/public/docs/
- [ ] PDF metadata embedded (Title, Author, Keywords, DOI)
- [ ] PDF URL references point to `/public/docs/` path

### Citations & Metadata
- [ ] All citation files have correct DOI: 10.5281/zenodo.17982540
- [ ] All citation files have correct author: Matthew Lukin Smawfield
- [ ] All citation files have correct title
- [ ] Keywords consistent across CITATION.cff, codemeta.json, index.html

### Git & Release
- [ ] `version.txt` updated with 0-6 high-level bullets (latest changes only, focus on pipeline/manuscript)
- [ ] Commit message prepared (see Phase 8 output)
- [ ] All changes staged/committed

## Files That Must Be Checked

### Core Metadata
- `VERSION.json`
- `version.txt`
- `README.md`
- `zenodo.txt` (must include: abstract, Website, Code Availability, Keywords)

### Site Files
- `site/index.html` (meta tags, PDF links, header)
- `site/manifest.json`
- `site/package.json`
- `site/CITATION.cff`
- `site/CITATION.bib` (check version in note field)
- `site/citation.json` (check version in note field)
- `site/codemeta.json`
- `site/components/abstract.html`

### Generated/Distribution
- `site/dist/index.html` (rebuild via `node build.js`)
- `13-TEP-WB-vX.Y-Codename.pdf` (root)
- `site/public/docs/13-TEP-WB-vX.Y-Codename.pdf`
- `13-TEP-WB-vX.Y-Codename.md` (auto-generated from site - will be rebuilt)

## Common Issues & Fixes

### Version String Not Updating
```bash
# Find all occurrences
grep -r "v0\.[0-9]" --include="*.html" --include="*.json" --include="*.md" .

# Replace carefully (example for v0.4 -> v0.5)
sed -i '' 's/v0\.4/v0.5/g' file.txt
```

### Abstract Mismatch
Always copy from `site/components/abstract.html` to all other locations. This is the source of truth.

### Title Mismatch
The title in `site/index.html` (inside `<title>` tags) is the source of truth. Update all other files to match:

```bash
# Extract the canonical title from index.html
TITLE="White Dwarf Cooling and the Temporal Equivalence Principle"

# Update files (example)
sed -i '' "s/^title: \".*\"/title: \"$TITLE\"/" site/CITATION.cff
sed -i '' "s/\"title\": \".*\"/\"title\": \"$TITLE\"/" site/citation.json
sed -i '' "s/\"name\": \".*\"/\"name\": \"$TITLE\"/" site/codemeta.json
```

**Important**: Title must match exactly including capitalization, spacing, and punctuation across all files.

### zenodo.txt Missing Template Sections
The zenodo.txt must include the full template with Website, Code Availability, Keywords, and Open Science Statement:

```bash
# Example zenodo.txt for TEP-WB:
cat > zenodo.txt << 'EOF'
{Abstract text here...}

Website: https://mlsmawfield.com/tep/wb
Code Availability: https://github.com/matthewsmawfield/TEP-WB

Keywords: gravitational lensing – dark matter – modified gravity – cosmology: theory – galaxies: kinematics and dynamics – temporal equivalence principle

Open Science Statement: This work is a preprint and is open to community review, ideas, and collaboration. All materials required for full reproducibility—including data downloads, analysis scripts, code, and manuscripts—are open-source. Feedback and contributions to further test these results are welcome.
EOF
```

### CITATION.bib or citation.json Version Outdated
These files often get missed during version bumps:

```bash
# Check current versions
grep "note.*Version" site/CITATION.bib
grep '"note":' site/citation.json

# Update to current version (e.g., 0.5)
sed -i '' 's/Version 0\.[0-9]/Version 0.5/g' site/CITATION.bib
sed -i '' 's/Version [0-9]\.[0-9]/Version 0.5/g' site/citation.json
```

### Old PDF Versions Still Present
Remove old PDF versions after generating new one:

```bash
# List all PDFs
ls -la *.pdf site/public/docs/*.pdf

# Remove old versions (keep only current)
rm "13-TEP-WB-v0.4-Tortola.pdf" 2>/dev/null
rm "site/public/docs/13-TEP-WB-v0.4-Tortola.pdf" 2>/dev/null
```

### PDF Not Found After Generation
- Check that `node build.js` ran successfully
- Verify `site/dist/index.html` exists
- Check Python script output for errors
- Verify PDF filename uses dash: `vX.Y-Codename.pdf`

### Date Confusion
Remember: **Citation dates = First published (Dec 2025)** | **Last updated = Today**

## Post-Publication Steps

After running this workflow and committing:

1. Create Git tag for the release:
   ```bash
   git tag -a vX.Y -m "vX.Y (Codename) - Release"
   git push origin vX.Y
   ```

2. Upload to Zenodo (if new DOI needed):
   - Use the generated PDF from root directory
   - Paste abstract from `zenodo.txt`
   - Verify metadata matches CITATION.cff

3. Verify live site:
   - Check https://matthewsmawfield.github.io/TEP-WB/
   - Verify version string in header
   - Verify PDF download works
   - Verify meta tags in page source
