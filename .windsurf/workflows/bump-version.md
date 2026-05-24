---
description: Bump version number across all project files and regenerate PDF
---

# Bump Version Workflow

Bump the project version number (e.g., v0.2 → v0.3), update all metadata, and regenerate the PDF.

## CRITICAL: Date Rules

There are TWO dates that serve different purposes:

1. **First Published Date** (date-released / datePublished) — NEVER CHANGES after initial release
   - This is the date the manuscript was first made public
   - In CITATION.cff: `date-released: 2026-03-19` (v0.1 initial release)
   - In codemeta.json: `datePublished: 2026-03-19`
   - In manifest.json: `first_published: "19 March 2026"`

2. **Last Updated Date** (dateModified / date) — CHANGES with every version bump
   - This is today's date when making the new version
   - In VERSION.json: `date: 2026-04-27`
   - In manifest.json: `date: 2026-04-27`
   - In codemeta.json: `dateModified: 2026-04-27`

**Rule**: When bumping version, ONLY update the "last updated" date. The "first published" date stays fixed forever.

## Overview

This workflow performs a complete version bump across the entire codebase:
1. Updates version number in all metadata files
2. Updates "last updated" date to TODAY
3. Leaves "first published" date unchanged
4. Regenerates the manuscript markdown
5. Copies to manuscripts/ folder

## Prerequisites

- Ensure all changes are committed or staged
- Verify the current version in `CITATION.cff`
- Know the new version number (e.g., 0.3)
- Codename stays the same (Kilifi) unless explicitly changing

## Quick Start

```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB"

# Build the site and regenerate markdown
npm run build:markdown --prefix site

# Copy to manuscripts folder
cp 13-TEP-WB-v0.3-Kilifi.md manuscripts/
```

## Manual Steps

### 1. Check Current Version

```bash
grep "^version:" CITATION.cff
grep "^date-released:" CITATION.cff
grep '"version"' VERSION.json
```

### 2. Update Version in All Files

Update these files with the NEW version (e.g., v0.2 → v0.3):

**Root level:**
- `CITATION.cff` — `version: "v0.3 (Kilifi)"` (KEEP date-released as first published)
- `VERSION.json` — `"version": "v0.3 (Kilifi)"`, update `"date": "2026-04-27"` (today)
- `README.md` — Update version string

**Site level:**
- `site/manifest.json` — `"version": "v0.3 (Kilifi)"`, update `"date": "2026-04-27"`
- `site/codemeta.json` — `"version": "v0.3 (Kilifi)"`, update `"dateModified": "2026-04-27"` (KEEP datePublished)
- `site/CITATION.cff` — `version: "v0.3 (Kilifi)"`
- `site/citation.json` — Update version

**Dist (auto-copied but check):**
- `site/dist/manifest.json`
- `site/dist/codemeta.json`
- `site/dist/CITATION.cff`
- `site/dist/citation.json`

### 3. Date Update Summary Table

| File | Field | Value | Action |
|------|-------|-------|--------|
| CITATION.cff | `version` | v0.3 (Kilifi) | UPDATE |
| CITATION.cff | `date-released` | 2026-03-19 | KEEP (first published) |
| VERSION.json | `version` | v0.3 (Kilifi) | UPDATE |
| VERSION.json | `date` | 2026-04-27 | UPDATE to today |
| site/manifest.json | `version` | v0.3 (Kilifi) | UPDATE |
| site/manifest.json | `date` | 2026-04-27 | UPDATE to today |
| site/manifest.json | `first_published` | "19 March 2026" | KEEP |
| site/codemeta.json | `version` | v0.3 (Kilifi) | UPDATE |
| site/codemeta.json | `datePublished` | 2026-03-19 | KEEP (first published) |
| site/codemeta.json | `dateModified` | 2026-04-27 | UPDATE to today |

### 4. Build and Generate

```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB"

# Build the site (updates dist/ and generates markdown)
npm run build:markdown --prefix site

# Verify the markdown header is correct
head -5 13-TEP-WB-v0.3-Kilifi.md

# Copy to manuscripts folder
cp 13-TEP-WB-v0.3-Kilifi.md manuscripts/
```

### 5. Verify Version Updates

Check all files were updated correctly:
```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-WB"

# Check version strings
grep -r "v0\.3" CITATION.cff VERSION.json site/manifest.json site/codemeta.json

# Check dates
grep "date-released" CITATION.cff
grep '"date"' VERSION.json site/manifest.json
grep "first_published" site/manifest.json
grep "datePublished\|dateModified" site/codemeta.json

# Verify markdown was generated with correct version
head -5 13-TEP-WB-v0.3-Kilifi.md
```

### 6. Git Status Check

```bash
git status
git diff --stat
```

Expected changed files:
- `CITATION.cff`
- `VERSION.json`
- `README.md`
- `site/manifest.json`
- `site/codemeta.json`
- `site/CITATION.cff`
- `site/citation.json`
- `site/dist/*` (auto-generated)
- `13-TEP-WB-v0.3-Kilifi.md` (auto-generated)
- `manuscripts/13-TEP-WB-v0.3-Kilifi.md` (copied)

## Commit Message Format

```
v0.3 (Kilifi) - Brief summary of changes

- Change 1
- Change 2
- Change 3
```

Or simply:
```
v0.3 (Kilifi) - Date and metadata updates

- Updated last modified date to 27 April 2026
- Version bump from v0.2 to v0.3
```

## Post-Bump Checklist

- [ ] Version updated in ALL metadata files
- [ ] **First published date UNCHANGED** (2026-03-19)
- [ ] **Last updated date set to TODAY**
- [ ] Manuscript markdown regenerated with correct header
- [ ] Markdown copied to manuscripts/ folder
- [ ] Git diff shows expected files changed
- [ ] Commit message prepared

## Common Mistakes to Avoid

1. **DON'T change `date-released` in CITATION.cff** — this is the FIRST publication date
2. **DON'T change `datePublished` in codemeta.json** — this is also the first publication date
3. **DO update `date` in VERSION.json** — this is the last modified date
4. **DO update `dateModified` in codemeta.json** — this is the last modified date
5. **DO update `date` in manifest.json** — this is the last modified date
